# -*- coding: utf-8 -*-
"""公文/学术稿件通用样式（方正小标宋标题、仿宋正文、三线表、页脚页码）"""
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement


def new_doc(body_font="仿宋_GB2312", body_size=16, line=1.5):
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.top_margin, sec.bottom_margin = Cm(3.0), Cm(2.5)
    sec.left_margin, sec.right_margin = Cm(3.0), Cm(2.5)
    st = doc.styles["Normal"]
    st.font.name = body_font
    st.font.size = Pt(body_size)
    st._element.rPr.rFonts.set(qn("w:eastAsia"), body_font)
    st._element.rPr.rFonts.set(qn("w:ascii"), "Times New Roman")
    st._element.rPr.rFonts.set(qn("w:hAnsi"), "Times New Roman")
    st.paragraph_format.line_spacing = line
    st.paragraph_format.space_before = Pt(0)
    st.paragraph_format.space_after = Pt(0)
    return doc


def set_run(r, font="仿宋_GB2312", size=16, bold=False):
    r.font.name = font
    r.font.size = Pt(size)
    r.bold = bold
    r._element.rPr.rFonts.set(qn("w:eastAsia"), font)
    if font not in ("宋体", "仿宋_GB2312", "楷体_GB2312", "黑体", "方正小标宋简体"):
        r._element.rPr.rFonts.set(qn("w:ascii"), "Times New Roman")
        r._element.rPr.rFonts.set(qn("w:hAnsi"), "Times New Roman")


class Doc:
    def __init__(self, body_font="仿宋_GB2312", body_size=16, line=1.5):
        self.doc = new_doc(body_font, body_size, line)
        self.font = body_font
        self.size = body_size
        self.line = line

    def title(self, text, size=22):
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(10)
        p.paragraph_format.line_spacing = 1.4
        set_run(p.add_run(text), "方正小标宋简体", size)
        return p

    def sub(self, text, size=14):
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_after = Pt(14)
        set_run(p.add_run(text), "楷体_GB2312", size)
        return p

    def h1(self, text):
        p = self.doc.add_paragraph()
        p.paragraph_format.line_spacing = self.line
        p.paragraph_format.space_before = Pt(12)
        p.paragraph_format.space_after = Pt(5)
        p.paragraph_format.keep_with_next = True
        set_run(p.add_run(text), "黑体", self.size)
        return p

    def h2(self, text):
        p = self.doc.add_paragraph()
        p.paragraph_format.line_spacing = self.line
        p.paragraph_format.space_before = Pt(7)
        p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.keep_with_next = True
        set_run(p.add_run(text), "楷体_GB2312", self.size, True)
        return p

    def p(self, text="", indent=2, size=None, font=None, bold=False,
          align=WD_ALIGN_PARAGRAPH.JUSTIFY, before=0, after=0):
        size = size or self.size
        p = self.doc.add_paragraph()
        p.alignment = align
        pf = p.paragraph_format
        pf.line_spacing = self.line
        pf.space_before = Pt(before)
        pf.space_after = Pt(after)
        if indent:
            pf.first_line_indent = Pt(size * indent)
        set_run(p.add_run(text), font or self.font, size, bold)
        return p

    def bullet(self, text, size=None):
        size = size or self.size
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        pf = p.paragraph_format
        pf.line_spacing = self.line
        pf.left_indent = Pt(size * 2)
        pf.first_line_indent = Pt(-size * 1.0)
        set_run(p.add_run(text), self.font, size)
        return p

    def tbl_caption(self, text):
        p = self.doc.add_paragraph()
        p.paragraph_format.line_spacing = 1.2
        p.paragraph_format.space_before = Pt(9)
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.keep_with_next = True
        set_run(p.add_run(text), "黑体", 10.5)

    def fig_caption(self, text):
        p = self.doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.line_spacing = 1.2
        p.paragraph_format.space_before = Pt(4)
        p.paragraph_format.space_after = Pt(10)
        set_run(p.add_run(text), "黑体", 10.5)

    def note(self, text, size=10):
        p = self.doc.add_paragraph()
        p.paragraph_format.line_spacing = 1.15
        p.paragraph_format.space_before = Pt(2)
        p.paragraph_format.space_after = Pt(8)
        set_run(p.add_run(text), "宋体", size)

    def table(self, headers, rows, widths=None, font_size=9.5, center=True):
        t = self.doc.add_table(rows=1, cols=len(headers))
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        t.autofit = False
        for j, h in enumerate(headers):
            c = t.rows[0].cells[j]; c.text = ""
            p = c.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.add_run(h)
        for r in rows:
            cells = t.add_row().cells
            for j, v in enumerate(r):
                cells[j].text = ""
                p = cells[j].paragraphs[0]
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER if (center and j > 0) else WD_ALIGN_PARAGRAPH.LEFT
                p.add_run("" if v is None else str(v))
        if widths:
            total = sum(widths); avail = 15.4
            if total > avail:
                widths = [w * avail / total for w in widths]
            lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed")
            t._tbl.tblPr.append(lay)
            for j, w in enumerate(widths):
                t.columns[j].width = Cm(w)
                for row in t.rows:
                    row.cells[j].width = Cm(w)
        trPr = t.rows[0]._tr.get_or_add_trPr()
        hdr = OxmlElement("w:tblHeader"); hdr.set(qn("w:val"), "true")
        trPr.append(hdr)
        self._three_line(t, font_size)
        # 短表整体连排、行内不断裂；长表允许跨页（表头已设重复）
        if len(t.rows) <= 16:
            for i, row in enumerate(t.rows):
                trPr = row._tr.get_or_add_trPr()
                cs = OxmlElement("w:cantSplit"); trPr.append(cs)
                if i < len(t.rows) - 1:
                    for cell in row.cells:
                        for p in cell.paragraphs:
                            p.paragraph_format.keep_with_next = True
        return t

    def _border(self, cell, edge, sz=8, val="single"):
        tcPr = cell._tc.get_or_add_tcPr()
        b = tcPr.find(qn("w:tcBorders"))
        if b is None:
            b = OxmlElement("w:tcBorders"); tcPr.append(b)
        el = b.find(qn("w:" + edge))
        if el is None:
            el = OxmlElement("w:" + edge); b.append(el)
        el.set(qn("w:val"), val); el.set(qn("w:sz"), str(sz)); el.set(qn("w:color"), "000000")

    def _three_line(self, table, font_size):
        for i, row in enumerate(table.rows):
            for cell in row.cells:
                for e in ("top", "bottom", "left", "right"):
                    self._border(cell, e, 0, "none")
                cell.vertical_alignment = 1
                for p in cell.paragraphs:
                    p.paragraph_format.line_spacing = 1.0
                    p.paragraph_format.space_before = Pt(1.5)
                    p.paragraph_format.space_after = Pt(1.5)
                    for r in p.runs:
                        set_run(r, "宋体", font_size, bold=(i == 0))
        for cell in table.rows[0].cells:
            self._border(cell, "top", 12); self._border(cell, "bottom", 6)
        for cell in table.rows[-1].cells:
            self._border(cell, "bottom", 12)

    def pagebreak(self):
        self.doc.add_page_break()

    def footer_pagenum(self):
        p = self.doc.sections[0].footer.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run()
        fld = OxmlElement("w:fldSimple")
        fld.set(qn("w:instr"), "PAGE")
        r._element.addnext(fld)
        for run in p.runs:
            set_run(run, "宋体", 10.5)

    def save(self, path):
        self.doc.save(str(path))
        return path
