"""Minimal streaming reader for STRICT-conformance OOXML (openpyxl refuses these)."""
import zipfile, re
from xml.etree import ElementTree as ET

def _sst(z):
    try: raw = z.read("xl/sharedStrings.xml")
    except KeyError: return []
    out=[]
    for si in ET.fromstring(raw):
        out.append("".join(t.text or "" for t in si.iter() if t.tag.endswith("}t")))
    return out

def sheet_map(z):
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    R = {}
    for rel in ET.fromstring(z.read("xl/_rels/workbook.xml.rels")):
        R[rel.get("Id")] = rel.get("Target")
    out={}
    for sh in wb.iter():
        if sh.tag.endswith("}sheet"):
            rid = [v for k,v in sh.attrib.items() if k.endswith("}id")][0]
            out[sh.get("name")] = "xl/" + R[rid].lstrip("/").replace("xl/","",1)
    return out

def rows(path, sheet, maxrows=None):
    z = zipfile.ZipFile(path); sst = _sst(z); target = sheet_map(z)[sheet]
    with z.open(target) as fh:
        cur, buf = None, {}
        for ev, el in ET.iterparse(fh, events=("end",)):
            if el.tag.endswith("}c"):
                ref = el.get("r") or ""
                col = re.match(r"([A-Z]+)", ref)
                ci = 0
                for ch in (col.group(1) if col else "A"): ci = ci*26 + ord(ch)-64
                t = el.get("t")
                v = el.find("{*}v"); isel = el.find("{*}is")
                if t == "inlineStr" or isel is not None:
                    val = "".join(x.text or "" for x in (isel.iter() if isel is not None else []) if x.tag.endswith("}t"))
                elif v is None: val = None
                elif t == "s": val = sst[int(v.text)]
                else: val = v.text
                buf[ci] = val
                el.clear()
            elif el.tag.endswith("}row"):
                n = max(buf) if buf else 0
                yield [buf.get(i) for i in range(1, n+1)]
                buf = {}; el.clear()
                if maxrows is not None:
                    maxrows -= 1
                    if maxrows <= 0: return
