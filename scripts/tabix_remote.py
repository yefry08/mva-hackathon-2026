#!/usr/bin/env python3
"""Cliente tabix por HTTP: consulta regiones de un VCF bgzip remoto sin bajarlo.

gnomAD v4.1 publica cada cromosoma como un .vcf.bgz de varios GB junto a un
indice .tbi de unos 50 KB, y el servidor acepta rangos de bytes. Con eso se puede
traer solo los bloques BGZF que cubren los genes que nos interesan: unos pocos MB
en lugar de 7.4 GB por cromosoma.

Esto tambien resuelve un problema de datos, no solo de ancho de banda: la consulta
que sale de esta maquina es un rango de bytes derivado de las coordenadas de un
gen del panel, no una posicion del paciente. Nada especifico del nino viaja.

Formato del indice: el mismo esquema de bins del BAM. Cabecera TBI\\1, por cada
secuencia una lista de bins con sus chunks (offsets virtuales de 64 bits, 48 bits
de offset de bloque comprimido y 16 de offset dentro del bloque descomprimido) y
un indice lineal con el offset minimo por ventana de 16 kb.

Uso como libreria:
    from tabix_remote import RemoteTabix
    tbx = RemoteTabix(url)
    for line in tbx.query("chr15", 40_160_000, 40_170_000):
        ...
"""
from __future__ import annotations

import gzip
import io
import struct
import sys
import urllib.request
import zlib

MAX_BIN = ((1 << 18) - 1) // 7
TAD_LIDX_SHIFT = 14


def reg2bins(beg: int, end: int):
    """Bins del esquema UCSC/BAM que pueden solapar [beg, end)."""
    end -= 1
    yield 0
    for start, shift in ((1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)):
        for k in range(start + (beg >> shift), start + (end >> shift) + 1):
            yield k


class RemoteTabix:
    def __init__(self, url: str, tbi_url: str | None = None, timeout: int = 60):
        self.url = url
        self.tbi_url = tbi_url or (url + ".tbi")
        self.timeout = timeout
        self._index = None

    # -- HTTP ---------------------------------------------------------------
    def _fetch(self, start: int, end: int) -> bytes:
        req = urllib.request.Request(self.url,
                                     headers={"Range": "bytes=%d-%d" % (start, end)})
        with urllib.request.urlopen(req, timeout=self.timeout) as r:
            return r.read()

    # -- indice -------------------------------------------------------------
    def _load_index(self):
        if self._index is not None:
            return self._index
        with urllib.request.urlopen(self.tbi_url, timeout=self.timeout) as r:
            raw = gzip.decompress(r.read())
        buf = io.BytesIO(raw)

        magic = buf.read(4)
        if magic != b"TBI\x01":
            raise ValueError("no es un indice tabix: %r" % magic)
        n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = \
            struct.unpack("<8i", buf.read(32))
        names = buf.read(l_nm).split(b"\x00")
        names = [n.decode() for n in names if n]

        refs = []
        for _ in range(n_ref):
            n_bin, = struct.unpack("<i", buf.read(4))
            bins = {}
            for _ in range(n_bin):
                bin_id, n_chunk = struct.unpack("<Ii", buf.read(8))
                chunks = []
                for _ in range(n_chunk):
                    cbeg, cend = struct.unpack("<QQ", buf.read(16))
                    chunks.append((cbeg, cend))
                bins[bin_id] = chunks
            n_intv, = struct.unpack("<i", buf.read(4))
            intervals = list(struct.unpack("<%dQ" % n_intv, buf.read(8 * n_intv))) \
                if n_intv else []
            refs.append({"bins": bins, "intervals": intervals})

        self._index = {"names": names, "refs": refs, "col_beg": col_beg}
        return self._index

    # -- consulta -----------------------------------------------------------
    def query(self, chrom: str, start: int, end: int):
        """start y end en coordenadas 1-based inclusivas, como en un VCF."""
        idx = self._load_index()
        names = idx["names"]
        if chrom not in names:
            alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            if alt not in names:
                raise KeyError("el indice no tiene %r (tiene %s...)" % (chrom, names[:3]))
            chrom = alt
        ref = idx["refs"][names.index(chrom)]

        beg0 = max(0, start - 1)
        chunks = []
        for b in reg2bins(beg0, end):
            chunks.extend(ref["bins"].get(b, ()))
        if not chunks:
            return

        # El indice lineal descarta chunks que terminan antes del inicio de la region.
        ivs = ref["intervals"]
        min_off = ivs[beg0 >> TAD_LIDX_SHIFT] if (beg0 >> TAD_LIDX_SHIFT) < len(ivs) else 0
        chunks = [(cb, ce) for cb, ce in chunks if ce > min_off]
        if not chunks:
            return

        # Unimos los rangos de bloques comprimidos para pedir pocos GET grandes.
        spans = sorted(((cb >> 16, ce >> 16) for cb, ce in chunks))
        merged = []
        for cb, ce in spans:
            if merged and cb <= merged[-1][1] + 65536:
                merged[-1][1] = max(merged[-1][1], ce)
            else:
                merged.append([cb, ce])

        for cb, ce in merged:
            data = self._fetch(cb, ce + 65535)
            for line in self._inflate(data):
                if not line or line.startswith("#"):
                    continue
                f = line.split("\t", 2)
                if len(f) < 2:
                    continue
                if f[0] != chrom and f[0] != chrom.lstrip("chr"):
                    continue
                try:
                    pos = int(f[1])
                except ValueError:
                    continue
                if start <= pos <= end:
                    yield line

    @staticmethod
    def _inflate(data: bytes):
        """Descomprime bloques BGZF concatenados, tolerando el ultimo truncado."""
        out = bytearray()
        while data:
            try:
                d = zlib.decompressobj(31)
                out.extend(d.decompress(data))
                if d.unused_data == data:
                    break
                data = d.unused_data
            except zlib.error:
                break
        text = out.decode("utf-8", errors="replace")
        # La primera y la ultima linea pueden venir cortadas por el rango pedido.
        lines = text.split("\n")
        return lines[1:-1] if len(lines) > 2 else []


def main() -> int:
    if len(sys.argv) != 5:
        print(__doc__)
        print("prueba: python scripts/tabix_remote.py URL chrom inicio fin")
        return 1
    url, chrom, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    n = 0
    for line in RemoteTabix(url).query(chrom, start, end):
        n += 1
        if n <= 3:
            print(line[:160])
    print("... %d registros en %s:%d-%d" % (n, chrom, start, end))
    return 0


if __name__ == "__main__":
    sys.exit(main())
