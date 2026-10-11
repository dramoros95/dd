# Builds the PDF discussion of the BBK3CW-Z NH3 results vs El Mir et al. (Ee = 74 eV)
import os
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image, Table,
                                TableStyle, PageBreak, KeepTogether)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

F = '/usr/share/fonts/truetype/dejavu/'
pdfmetrics.registerFont(TTFont('DV', F + 'DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DVB', F + 'DejaVuSans-Bold.ttf'))
pdfmetrics.registerFontFamily('DV', normal='DV', bold='DVB', italic='DV', boldItalic='DVB')

D = os.path.dirname(os.path.abspath(__file__)) + '/../'
OUT = D + 'report/NH3_BBK3CWZ_vs_ElMir_Ee74eV_discussion.pdf'

ss = getSampleStyleSheet()
body = ParagraphStyle('b', parent=ss['Normal'], fontName='DV', fontSize=9.6, leading=13.4, spaceAfter=5)
bul = ParagraphStyle('u', parent=body, leftIndent=12, bulletIndent=2, spaceAfter=3)
h1 = ParagraphStyle('h1', parent=body, fontName='DVB', fontSize=16, leading=20, spaceAfter=4)
h2 = ParagraphStyle('h2', parent=body, fontName='DVB', fontSize=12, leading=16, spaceBefore=8, spaceAfter=4)
sub = ParagraphStyle('s', parent=body, fontSize=9, textColor=colors.HexColor('#52514e'))
cap = ParagraphStyle('c', parent=body, fontSize=8.6, leading=11.5, textColor=colors.HexColor('#3a3a38'))
cell = ParagraphStyle('t', parent=body, fontSize=8.4, leading=10.5, spaceAfter=0)
cellb = ParagraphStyle('tb', parent=cell, fontName='DVB')

P = lambda t, s=body: Paragraph(t, s)
B = lambda t: Paragraph(t, bul, bulletText='•')

def table(rows, widths, header=True):
    data = [[Paragraph(c, cellb if (header and i == 0) else cell) for c in r] for i, r in enumerate(rows)]
    t = Table(data, colWidths=widths, repeatRows=1 if header else 0)
    st = [('GRID', (0, 0), (-1, -1), 0.4, colors.HexColor('#bdbcb8')),
          ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
          ('TOPPADDING', (0, 0), (-1, -1), 2.5), ('BOTTOMPADDING', (0, 0), (-1, -1), 2.5)]
    if header:
        st.append(('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#ecebe7')))
    t.setStyle(TableStyle(st))
    return t

def fig(name, w=12.3 * cm):
    im = Image(D + name)
    im.drawHeight = w * im.imageHeight / im.imageWidth
    im.drawWidth = w; im.hAlign = "CENTER"
    return im

a1, a3, e1 = '3a<sub>1</sub>', '2a<sub>1</sub>', '1e'
s = []
s += [P('BBK3CW-Z results for NH<sub>3</sub> vs experiment', h1),
      P('El Mir et al., J. Phys. B 54, 015201 (2020) — E<sub>s</sub> = 500 eV, θ<sub>s</sub> = −6°, '
        'E<sub>e</sub> = 74 eV, coplanar asymmetric geometry. Orbitals ' + a1 + ', 1e, ' + a3 +
        ' (Moccia one-centre SCF wave functions).', sub), Spacer(1, 6)]

s += [P('How to read the figures', h2),
      B('Every model curve is scaled to the maximum of the experimental binary peak. The comparison '
        'is therefore of <b>shapes and binary/recoil ratios</b>, not of absolute magnitudes.'),
      B('The experiment has no points between 160° and 200°, the angular range the detectors cannot '
        'reach. The experimental recoil maximum is taken in 200–335°.'),
      B('The momentum transfer <b>K</b> points at θ<sub>K</sub> ≈ 50° (46° for ' + a3 + '), so −<b>K</b> '
        'is at ≈ 230°. In a first-order picture the binary peak lies near <b>K</b> and the recoil peak '
        'near −<b>K</b>.'),
      B('Three versions of the model are shown: the original BBK3CW-Z (Fig. 1); with the '
        'short-range projectile–core potential, option isr = 1 (Fig. 2); and with the ejected wave '
        'orthogonalized to the orbital, option iorth = 1, with and without isr (Fig. 3).')]

s += [P('Summary of the numbers', h2),
      table([['Orbital / version', 'Binary peak', 'Recoil peak', 'Binary / recoil'],
             ['<b>' + a1 + '</b> experiment', '54°', '242°', '2.1'],
             ['BBK3CW-Z', '60°', 'flat (max 198°)', '14.8'],
             ['+ short-range (isr = 1)', '66°', '252°', '0.98'],
             ['orthogonalized', '54°', 'none', '30'],
             ['<b>1e</b> experiment', '49°', '240°', '2.5'],
             ['BBK3CW-Z', '54°', 'flat (max 192°)', '23.4'],
             ['+ short-range', '60°', 'broad plateau', '1.8'],
             ['orthogonalized', '54°', 'none', '37'],
             ['<b>' + a3 + '</b> experiment', '50°', '230°', '1.1'],
             ['BBK3CW-Z', '48°', '246° (weak)', '5.3'],
             ['+ short-range', '66° (deformed)', '258°', '0.7'],
             ['orthogonalized', '42°', 'none', '3.0']],
            [5.4 * cm, 3.2 * cm, 3.8 * cm, 3.2 * cm]),
      Spacer(1, 3),
      P('“none”: for the orthogonalized runs the largest value between 180° and 340° is at 336°, on '
        'the tail of the binary peak; there is no real recoil peak.', cap)]

# ---- Figure 1
s += [PageBreak(), P('1. Original BBK3CW-Z vs experiment', h2),
      fig('compare_NH3_Ee74eV_vs_ElMir.png', 11.2 * cm),
      P('<b>Figure 1.</b> BBK3CW-Z (blue line) and the El Mir et al. experiment (orange diamonds) for '
        + a3 + ', 1e and ' + a1 + '. Each model curve is scaled to the experimental binary maximum.', cap),
      P('<b>What works</b>'),
      B('The <b>binary peak position</b> is close to experiment: 48–60° against 49–54°. The model is a few '
        'degrees too far forward for ' + a1 + ' and 1e.'),
      B('The <b>front side of the binary peak</b> (25–60°) follows the data well.'),
      P('<b>What fails</b>'),
      B('<b>The recoil peak is almost missing.</b> For ' + a1 + ' and 1e the model is a flat plateau 5–10 times '
        'below the data, with no maximum near −<b>K</b>. For ' + a3 + ' there is a weak bump at 246°, but it is '
        'still 4–5 times too low.'),
      B('<b>The binary peak is too narrow on its large-angle side.</b> At θ<sub>e</sub> = 120°, the 1e experiment is '
        'about 1×10<super>−16</super> and the model about 1.5×10<super>−17</super>, a factor of about 7; the same '
        'holds for ' + a1 + '. In the experiment the cross section falls slowly from the binary peak towards the '
        'recoil region; in the model it drops sharply after 90°.'),
      B('<b>' + a3 + ' has a dip at about 80° and a second bump at 105°</b> that the experiment does not show '
        '(it has a smooth shoulder). ' + a3 + ' is mostly N 2s, whose momentum distribution has a node, and the '
        'model shows that node too sharply. In the real collision the distortion of the ejected electron '
        'would smear it out.'),
      P('The first two failures have the same cause: in the model, the ejected electron is not scattered '
        'enough by the ion. Whatever leaves the binary direction (slow fall-off, valley filling, recoil '
        'peak) requires the ion to deflect the electron, and the BBK3CW-Z ejected wave barely does.')]

# ---- Figure 2
s += [PageBreak(), P('2. With the short-range term (isr = 1)', h2),
      fig('compare_NH3_Ee74eV_SR_vs_ElMir.png'),
      P('<b>Figure 2.</b> BBK3CW-Z (dashed blue), BBK3CW-Z with the short-range projectile–core potential '
        '−(Z<sub>ion</sub>(r)−1)/r (solid green), and experiment (orange diamonds).', cap),
      P('This option adds the projectile–core attraction inside the molecule.'),
      B('<b>1e</b> is the best result so far: the binary peak is right, and the curve from 270° to 335° '
        'falls on the data. But 130–200° is filled by a plateau that is too high, so the recoil looks like '
        'a broad shelf instead of a peak centred at 240°.'),
      B('<b>' + a1 + '</b>: the recoil is now about 2 times too high and much too broad; it starts right after '
        'the binary peak, at about 110°.'),
      B('<b>' + a3 + '</b>: the binary peak loses its shape and shifts to 66°. The curve is almost flat from 60° '
        'to 270°, which is unphysical. ' + a3 + ' has the largest s part, so the large overlap between the '
        'ejected Coulomb wave and the orbital (0.37) inflates this term the most.'),
      B('<b>Binary peak position</b>: it moves forward to 60–66° in all three orbitals, away from the '
        'experimental 49–54°.'),
      P('The term pushes in the right direction (it produces back-scattering) but in the wrong way: it '
        'adds a nearly isotropic contribution instead of a peak centred at −<b>K</b>. That isotropic part is '
        'what one obtains when the ejected wave is not orthogonal to the orbital.')]

# ---- Figure 3
s += [PageBreak(), P('3. Orthogonalized ejected wave (iorth = 1)', h2),
      fig('compare_NH3_Ee74eV_ORTH_vs_ElMir.png'),
      P('<b>Figure 3.</b> BBK3CW-Z (dashed blue), orthogonalized ejected wave (solid green), orthogonalized '
        'with the short-range term (dash-dot purple), and experiment (orange diamonds).', cap),
      B('<b>Binary peak</b>: unchanged and good (42–54°).'),
      B('<b>Recoil without the short-range term</b>: even lower than the original, with a deep minimum '
        'around 255–270°.'),
      B('<b>With the short-range term</b>: the isotropic shelf of Fig. 2 is gone, so ' + a3 + ' recovers its '
        'binary peak shape, but the recoil region falls back to about 10 times below the data.'),
      P('Most of the extra recoil produced by the short-range term was therefore the overlap artefact; the '
        'physical part left is far too weak.')]

# ---- Conclusion
s += [PageBreak(), P('What the three figures say together', h2),
      table([['Observation', 'Meaning'],
             ['Binary peak correct in all versions',
              'The initial state (Moccia orbitals), the first-order projectile treatment and the kinematics are fine'],
             ['Recoil and large-angle side 5–10 times low, for all 3 orbitals',
              'A missing physical mechanism, not an orbital-specific problem'],
             ['Short-range term raises it, but isotropically',
              'That increase comes mainly from the non-orthogonal overlap, not from real back-scattering'],
             ['Orthogonalization removes the increase',
              'Once the artefact is gone, almost no real back-scattering is left in the ejected wave']],
            [6.6 * cm, 9.0 * cm]),
      Spacer(1, 8),
      P('All four observations point to the <b>ejected-electron wave</b>. It does not contain the strong '
        'back-scattering of the NH<sub>3</sub><super>+</super> core. Solving Eq. (11) of the theory report '
        'exactly with the same Z(r) gives large short-range phase shifts (2.6, 1.3 and 0.5 rad for l = 0, 1, 2) '
        'and an elastic e<super>−</super> + NH<sub>3</sub><super>+</super> cross section 30–140 times the '
        'Rutherford (Z = 1) value at 120–180°. The variable-charge Coulomb wave of BBK3CW-Z keeps the Z = 1 '
        'phases and cannot reproduce this back-scattering. The proposed next step is therefore an exact '
        'distorted wave for the ejected electron (see theory_distorted_wave.md), not another adjustment of '
        'Z(r) or of the potential term.'),
      P('<b>Caveat.</b> This is a shape comparison only. If El Mir et al. put the three orbitals on a common '
        'relative scale, the orbital-to-orbital ratios would give an additional test.')]

def footer(c, d):
    c.saveState(); c.setFont('DV', 7.5); c.setFillColor(colors.HexColor('#8a8985'))
    c.drawRightString(A4[0] - 2 * cm, 1.2 * cm, 'page %d' % d.page)
    c.drawString(2 * cm, 1.2 * cm, 'NH3 BBK3CW-Z vs El Mir et al. (Ee = 74 eV)')
    c.restoreState()

SimpleDocTemplate(OUT, pagesize=A4, leftMargin=2 * cm, rightMargin=2 * cm, topMargin=1.6 * cm,
                  bottomMargin=1.8 * cm, title='NH3 BBK3CW-Z vs El Mir et al. (Ee = 74 eV)',
                  author='').build(s, onFirstPage=footer, onLaterPages=footer)
print(OUT)
