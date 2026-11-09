"""Build the submission DOCX (Word, editable equations) and PDF (via HTML + Chrome)."""
import re, subprocess, zipfile, shutil, os

NAME = "אמינות_רשתות_פרויקט_חלק_א"

# ---------- DOCX ----------
subprocess.run(["pandoc", "report.md", "-o", "tmp.docx", "--reference-doc=reference.docx"], check=True)
zin = zipfile.ZipFile("tmp.docx")
zout = zipfile.ZipFile(NAME + ".docx", "w", zipfile.ZIP_DEFLATED)
for it in zin.infolist():
    d = zin.read(it.filename)
    if it.filename == "word/document.xml":
        s = d.decode()

        def fix_par(m):
            p = m.group(0)
            if 'w:val="SourceCode"' in p or "<m:oMathPara" in p:
                p = p.replace("<w:bidi />", "")
                # left-align, LTR
                p = re.sub(r"(<w:pPr>)", r'\1<w:jc w:val="left" />', p, count=1)
            return p

        s = re.sub(r"<w:p>.*?</w:p>|<w:p .*?</w:p>", fix_par, s, flags=re.S)
        d = s.encode()
    if it.filename == "word/styles.xml":
        s = d.decode()
        # code style: monospace font
        s = s.replace('w:styleId="SourceCode"', 'w:styleId="SourceCode"')
        d = s.encode()
    zout.writestr(it, d)
zout.close()
os.remove("tmp.docx")

# ---------- PDF ----------
CSS = """
@page { size: A4; margin: 18mm 16mm; }
body { font-family: Arial, 'Liberation Sans', sans-serif; font-size: 11pt; line-height: 1.5; direction: rtl; color:#111; }
header#title-block-header { text-align: center; margin-bottom: 18px; border-bottom: 2px solid #1f3b63; padding-bottom: 10px;}
h1.title { font-size: 22pt; color:#1f3b63; margin: 0 0 6px 0;}
p.subtitle { font-size: 13pt; color:#1f3b63; margin: 0 0 6px 0;}
p.author, p.date { margin: 2px 0; color:#444;}
h1 { font-size: 15pt; color:#1f3b63; border-bottom: 1px solid #c9d3e0; padding-bottom: 3px; margin-top: 22px; page-break-after: avoid;}
h2 { font-size: 12.5pt; color:#1f3b63; margin-top: 14px; page-break-after: avoid;}
table { border-collapse: collapse; margin: 8px auto; page-break-inside: avoid; }
th, td { border: 1px solid #9aa8bb; padding: 3px 10px; text-align: center; direction: ltr; unicode-bidi: plaintext;}
th { background: #e8eef6; }
pre { direction: ltr; text-align: left; background:#f5f7fa; border:1px solid #dde3ea; padding: 8px 10px; font-size: 8.6pt; line-height: 1.35; page-break-inside: avoid;}
code { font-family: 'DejaVu Sans Mono', Consolas, monospace; }
math[display="block"] { direction: ltr; margin: 8px 0; font-size: 115%; }
figure { text-align:center; margin: 10px 0; page-break-inside: avoid;}
figure img { width: 46%; }
figcaption { font-size: 9.5pt; color:#555;}
li { margin: 2px 0; }
"""
open("style.css", "w").write(CSS)
subprocess.run(["pandoc", "report.md", "-s", "--mathml", "--embed-resources", "--css=style.css",
                "-o", "report.html", "--metadata", "pagetitle=report"], check=True)
subprocess.run(["google-chrome", "--headless=new", "--no-sandbox", "--disable-gpu",
                "--no-pdf-header-footer", f"--print-to-pdf={NAME}.pdf", "report.html"],
               check=True, capture_output=True)
print("built", NAME)
