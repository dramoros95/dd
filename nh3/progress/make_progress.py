# Builds the progress synthesis PDF: protocol objectives vs Q2/Q3 reports and the NH3 work
import os
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak, KeepTogether
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

F = '/usr/share/fonts/truetype/dejavu/'
pdfmetrics.registerFont(TTFont('DV', F + 'DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DVB', F + 'DejaVuSans-Bold.ttf'))
pdfmetrics.registerFontFamily('DV', normal='DV', bold='DVB', italic='DV', boldItalic='DVB')
H = os.path.dirname(os.path.abspath(__file__))
N = H + '/../'
OUT = H + '/Progress_3CWZ_CH4_NH3_protocol_status.pdf'

ss = getSampleStyleSheet()
body = ParagraphStyle('b', parent=ss['Normal'], fontName='DV', fontSize=9.4, leading=13, spaceAfter=4)
bul = ParagraphStyle('u', parent=body, leftIndent=12, bulletIndent=2, spaceAfter=2.5)
h1 = ParagraphStyle('h1', parent=body, fontName='DVB', fontSize=16, leading=20, spaceAfter=3)
h2 = ParagraphStyle('h2', parent=body, fontName='DVB', fontSize=12.5, leading=16, spaceBefore=9, spaceAfter=4)
h3 = ParagraphStyle('h3', parent=body, fontName='DVB', fontSize=10.5, leading=14, spaceBefore=6, spaceAfter=3)
sub = ParagraphStyle('s', parent=body, fontSize=9, textColor=colors.HexColor('#52514e'))
cap = ParagraphStyle('c', parent=body, fontSize=8.4, leading=11, textColor=colors.HexColor('#3a3a38'))
cell = ParagraphStyle('t', parent=body, fontSize=8.1, leading=10.2, spaceAfter=0)
cellb = ParagraphStyle('tb', parent=cell, fontName='DVB')
P = lambda t, s=body: Paragraph(t, s)
B = lambda t: Paragraph(t, bul, bulletText='•')
ST = {'DONE': '#d7efe3', 'PARTIAL': '#fbecc8', 'TO DO': '#f8d9d2'}

def table(rows, widths, status_col=None):
    data = [[Paragraph(c, cellb if i == 0 else cell) for c in r] for i, r in enumerate(rows)]
    t = Table(data, colWidths=widths, repeatRows=1)
    st = [('GRID', (0, 0), (-1, -1), 0.4, colors.HexColor('#bdbcb8')), ('VALIGN', (0, 0), (-1, -1), 'TOP'),
          ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#ecebe7')),
          ('TOPPADDING', (0, 0), (-1, -1), 2.5), ('BOTTOMPADDING', (0, 0), (-1, -1), 2.5)]
    if status_col is not None:
        for i, r in enumerate(rows[1:], 1):
            k = r[status_col].replace('<b>', '').replace('</b>', '')
            if k in ST: st.append(('BACKGROUND', (status_col, i), (status_col, i), colors.HexColor(ST[k])))
    t.setStyle(TableStyle(st))
    return t

def fig(path, w):
    im = Image(path); im.drawHeight = w * im.imageHeight / im.imageWidth; im.drawWidth = w; im.hAlign = 'CENTER'
    return im

a1, a3, t2 = '3a<sub>1</sub>', '2a<sub>1</sub>', '1t<sub>2</sub>'
NH3, CH4 = 'NH<sub>3</sub>', 'CH<sub>4</sub>'
s = []
s += [P('3CWZ for ' + CH4 + ' and ' + NH3 + ': progress against the protocol', h1),
      P('Synthesis of the research protocol (O. Zaidi, CRNA), the 2nd- and 3rd-quarter 2026 scientific reports, '
        'and the code and calculations of this work (repository dramoros95/dd, folder nh3/).', sub), Spacer(1, 4)]

s += [P('1. Sources', h2),
      table([['Source', 'Content', 'Role in this synthesis'],
             ['Protocol (Protocole_ZAIDI_VF)', 'Objectives, methodology, 12-month plan, expected results',
              'Reference checklist (sections 2–3)'],
             ['Q2 report (2nd quarter 2026)', 'Molecular 3CWZ formalism: kinematics, TDCS, 3CWZ final state, Moccia OCE, '
              'Euler average (Eq. 19), T<sub>dir</sub> = T<super>1</super> − T<super>0</super>, Fourier transform of C<sub>ab</sub>, '
              'split at r = d (Eqs. 23–36)', 'Theory of the amplitude'],
             ['Q3 report (3rd quarter 2026)', 'Z<sub>eff</sub>(r) from the monopole Hartree potential of the Moccia density; '
              'self-interaction correction for the ejected electron; Z(r) curves for ' + CH4 + ' and ' + NH3,
              'Theory of the variable charge'],
             ['This work', 'NH<sub>3</sub> code (nh3_opt_claude.for), runs vs El Mir et al. (2020), short-range and '
              'orthogonalization tests, ejected-wave diagnosis, tests under Rezkallah et al. (2011) conditions',
              'Implementation and results']],
            [3.6 * cm, 8.0 * cm, 4.6 * cm])]

# ---- objectives
s += [P('2. Protocol objectives: status', h2),
      table([['Protocol objective (§4)', 'Where it is covered', 'Status', 'What remains'],
             ['O1. Adapt 3CWZ to molecules with Moccia OCE orbitals; Z(r) for each MO',
              'Q2 (formalism); Q3 (Z(r) derivation); code: Moccia Table I for ' + NH3 + ' (' + a1 + ', 1e, ' + a3 + '), '
              'ZSCR = Q3 Eq. (16), analytic orientation average; ' + CH4 + ' ' + t2 + ' regression identical to the original code',
              'DONE', 'Per-orbital decomposition of Z(r) (Q3 §4.3), optional'],
             ['O2. TDCS for the orbitals of ' + CH4 + ' and ' + NH3 + ', several kinematics',
              NH3 + ': 3 orbitals at E<sub>s</sub> = 500 eV, θ<sub>s</sub> = −6°, E<sub>e</sub> = 37 and 74 eV; EMS at 1200 eV; '
              'DDCS at 500 eV. ' + CH4 + ': ' + t2 + ' with the original code',
              'PARTIAL', CH4 + ' ' + '2a<sub>1</sub> and new ' + CH4 + ' runs with the same code; other kinematics'],
             ['O3. PCI: compare results with and without C<sub>12</sub>', '—', 'TO DO',
              'Add a switch (C<sub>12</sub> = 1, N<sub>ab</sub> = 1) and rerun; cheap'],
             ['O4. Validate against experiment and theory: ' + CH4 + ' [35–37], ' + NH3 + ' [38–39]',
              NH3 + ': El Mir et al. 2020 (protocol ref. [9]) at 74 eV; also Bawagan EMS and Opal DDCS',
              'PARTIAL', 'Nixon et al. 2013 [39] (low energy, needs exchange); Rejoub et al. 2001 [38] (total CS, '
              'needs converged integrals); ' + CH4 + ' Lahmam-Bennani 2009 [35]'],
             ['Expected result: contribution of the variable charge',
              'Z(r) computed and cross-checked; its effect on the TDCS not isolated yet', 'PARTIAL',
              'Same runs with constant Z = 1 (3C/BBK limit) for a direct comparison'],
             ['Expected result: validity domain of 3CWZ for molecules',
              'Binary peak reproduced; recoil peak 5–10 times too small (' + NH3 + ', 74 eV); cause identified (ejected wave); '
              'numerical domain mapped (slow electrons, K &lt; 3 a.u.)', 'PARTIAL',
              'Extend to other kinematics; test the distorted-wave fix']],
            [3.4 * cm, 6.1 * cm, 1.95 * cm, 4.75 * cm], status_col=2),
      Spacer(1, 6),
      KeepTogether([P('Methodology steps (protocol §5)', h3),
      table([['Step (planned duration)', 'Status', 'Comment'],
             ['Bibliography (2 months)', 'PARTIAL', 'El Mir et al. 2020, Rezkallah et al. 2011, Moccia 1964 used; '
              'Champion et al. 2006 (1DW for H<sub>2</sub>O) and Nixon et al. 2013 still to read'],
             ['OCE orbitals and Z(r) in Fortran (1 month)', 'DONE', 'All MOs of Moccia Table I; norms = 1; Z(∞) = 1 (ion), 0 (neutral)'],
             ['Molecular 3CWZ code (1 month)', 'DONE', 'Plus options isr (short-range potential) and iorth (orthogonalization)'],
             ['Calculations ' + CH4 + ' with PCI study (3 months)', 'TO DO', 'Only the regression test was run here'],
             ['Calculations ' + NH3 + ' vs experiment (3 months)', 'PARTIAL', '74 eV complete; 37 eV computed, data not digitized'],
             ['Analysis and discussion (2 months)', 'PARTIAL', 'Recoil diagnosis and numerical limits (sections 5–6)']],
            [5.0 * cm, 1.95 * cm, 9.25 * cm], status_col=1)])]

# ---- what each report establishes + Z cross-check
s += [Spacer(1, 8), P('3. The variable charge: Q3 report vs code', h2),
      P('The Q3 note derives Z<sub>eff</sub>(r) as the monopole part of the Hartree potential of the Moccia density '
        '(its Eq. 16), with the ejected electron screened by one electron less (self-interaction correction, Eqs. 18–19). '
        'The code computes exactly this expression (function ZSCR). Values from the code:'),
      table([['r (a.u.)', '0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0', '4.0'],
             [NH3 + ' incident (neutral)', '7.000', '2.690', '1.113', '0.549', '0.541', '0.227', '0.091', '0.013'],
             [NH3 + ' ejected (' + a1 + ' hole)', '7.000', '3.094', '1.806', '1.402', '1.474', '1.198', '1.079', '1.011']],
            [3.8 * cm] + [1.55 * cm] * 8),
      Spacer(1, 4),
      fig(H + '/Zr_code_CH4_NH3.png', 15.5 * cm),
      P('<b>Figure 1.</b> Z<sub>eff</sub>(r) from the code (files Zr_*_code.dat).', cap),
      P('Every prediction of Q3 is confirmed by the code:'),
      B('Z(0) = Z<sub>N</sub> (6 for ' + CH4 + ', 7 for ' + NH3 + ').'),
      B('Asymptotes 0 (incident) and 1 (ejected) reached by r ≈ 4 a.u. (Z<sub>ion</sub>(4) = 1.011). This also justifies '
        'the cut-off d = 4 a.u. of Q2 Eq. (34).'),
      B('Z<sub>ej</sub>(r) ≥ Z<sub>inc</sub>(r) at every r.'),
      B('Local maximum exactly at r = a<sub>H</sub> (1.928 a.u. for ' + NH3 + ', 2.08 a.u. for ' + CH4 + '): Z = 1.53 / 0.61 '
        '(' + NH3 + ' ejected / incident). In the shell model it is a kink, from min(r/a<sub>H</sub>, 1).'),
      B('The three ' + NH3 + ' orbitals give nearly the same Z<sub>ion</sub>(r) (within 0.05).')]

# ---- results
s += [PageBreak(), P('4. Results for ' + NH3 + ' (this work)', h2),
      P('4.1 TDCS vs El Mir et al. (E<sub>s</sub> = 500 eV, θ<sub>s</sub> = −6°, E<sub>e</sub> = 74 eV)', h3),
      fig(N + 'elmir_conditions/compare_NH3_Ee74eV_vs_ElMir.png', 8.6 * cm),
      P('<b>Figure 2.</b> BBK3CW-Z (scaled to the experimental binary maximum) vs experiment. '
        'Detailed discussion: elmir_conditions/report/NH3_BBK3CWZ_vs_ElMir_Ee74eV_discussion.pdf.', cap),
      table([['Orbital', 'Exp. binary / recoil', 'BBK3CW-Z', '+ short-range (isr = 1)', 'orthogonalized'],
             [a1, '2.1', '14.8', '0.98', '30'], ['1e', '2.5', '23.4', '1.8', '37'], [a3, '1.1', '5.3', '0.7', '3.0']],
            [2.4 * cm, 3.4 * cm, 2.6 * cm, 4.2 * cm, 3.4 * cm]),
      Spacer(1, 4),
      B('<b>Binary peak</b>: position (48–60° vs 49–54°) and front side correct in every version.'),
      B('<b>Recoil peak and large-angle side</b>: 5–10 times too low for all three orbitals.'),
      B('<b>Short-range potential</b>: raises the recoil, but isotropically. Most of the rise comes from the '
        'non-orthogonality of the ejected wave; orthogonalization removes it and the recoil falls back.')]

s += [P('4.2 Diagnosis: the ejected-electron wave', h3),
      P('The ejected wave of Q2 Eq. (9) with the local α(r) = Z(r)/k does not solve Eq. (11) for the Q3 potential. '
        'The exact solution in the same Z<sub>ion</sub>(r) at 74 eV gives large short-range phase shifts, and the '
        'NH<sub>3</sub><super>+</super> core back-scatters far more than a Z = 1 field:'),
      table([['l', '0', '1', '2', '3', '4', '5'], ['δ<sub>l</sub> (rad)', '2.60', '1.34', '0.48', '0.26', '0.17', '0.09']],
            [3 * cm] + [1.6 * cm] * 6),
      Spacer(1, 3),
      table([['elastic angle', '30°', '60°', '90°', '120°', '150°', '180°'],
             ['exact / Rutherford (Z = 1)', '2.5', '4.8', '4.9', '29', '91', '137']], [4 * cm] + [1.5 * cm] * 6),
      Spacer(1, 3),
      P('A plane-wave estimate confirms that the recoil cannot come from the orbitals: ρ(q<sub>binary</sub>)/ρ(q<sub>recoil</sub>) '
        '= 38 (' + a1 + '), 59 (1e), 155 (' + a3 + ') against 1–2.5 measured. The recoil must come from the scattering of '
        'the ejected electron by the ion. Proposed fix: an exact distorted wave for the ejected electron in the Q3 '
        'potential (theory_distorted_wave.md).')]

s += [PageBreak(), P('4.3 Tests under the conditions of Rezkallah et al. (2011)', h3),
      fig(N + 'rezkallah_conditions/EMS_NH3_1200eV_symmetric_noncoplanar_45deg.png', 16 * cm),
      P('<b>Figure 3.</b> EMS momentum profiles (1200 eV, 45°, 600/600 eV) vs Bawagan et al. and the paper\'s PWBA.', cap),
      fig(N + 'rezkallah_conditions/DDCS_NH3_Ei500eV_thetae120deg.png', 9.5 * cm),
      P('<b>Figure 4.</b> DDCS at E<sub>i</sub> = 500 eV, θ<sub>e</sub> = 120° vs Opal et al.', cap),
      B('<b>EMS</b>: peak positions correct; shapes of ' + a3 + ' and ' + a1 + ' correct; the 1e high-momentum tail is '
        '~1.5 times low (orbital). Absolute value 0.57 of the PWBA with default grids (0.76 with finer grids for ' + a1 + ').'),
      B('<b>DDCS</b>: model/experiment = 1.79, 1.58, 1.35 at E<sub>e</sub> = 5, 20, 60 eV; 7.1 at 150 eV, where the excess is '
        'numerical (large momentum transfer; finer grids lower it by 12–30 times).'),
      B('<b>Numerical domain</b>: the default grids (nleg = 10, nd = 4, md = 1) are reliable to ~20% only for slow '
        'outgoing electrons and K &lt; 3 a.u. The El Mir conditions (K ≈ 0.85 a.u., 74 eV) are inside it.')]

# ---- consistency
s += [Spacer(1, 8), P('5. Points to correct or clarify in the documents', h2),
      table([['Document', 'Point', 'Suggested change'],
             ['Q3, table §3', NH3 + ' occupied MOs listed as "1A1, 2A1, 1A2, 1E"',
              'Moccia Table I: 1a<sub>1</sub>, 2a<sub>1</sub>, 3a<sub>1</sub>, 1e (no a<sub>2</sub> orbital); 3a<sub>1</sub> is the HOMO (lone pair)'],
             ['Q3, table §3 and §4.2', 'a<sub>H</sub>(N–H) = 1.91 a.u.',
              'Moccia Table I uses R = 1.928 a.u. (also in the code); use the value of the orbitals'],
             ['Q3, §3.1', '"this is how PCI enters"', 'The +1 tail is the electron–ion final-state interaction; PCI is the '
              'electron–electron factor C<sub>ab</sub> (Q2 Eq. 13)'],
             ['Q3, Figs. 1–2', 'Colours of incident/ejected swapped between the two captions', 'Use the same convention'],
             ['Q2, Eq. (2.25)/(22)', 'σ ∝ |T<sub>dir</sub>|² + |T<sub>exc</sub>|² + |T<sub>dir</sub> − T<sub>exc</sub>|²',
              'Check the factor ½ of the spin average: ½(...) = |T<sub>d</sub>|² + |T<sub>e</sub>|² − Re(T<sub>d</sub>T<sub>e</sub>*)'],
             ['Q2, Eq. (30)', 'Prefactor written 1/2π³', 'Should read 1/(2π)³ as in Eq. (29)'],
             ['Q2, Eq. (7) vs Q3', 'Projectile–core potential −1/r<sub>0</sub> (frozen core), while the waves use Z(r)',
              'State the choice; the full −Z<sub>ion</sub>(r)/r<sub>0</sub> is option isr = 1 (tested, section 4.1)'],
             ['Q2, §3.4', 'Exchange amplitude dropped', 'Fine for E<sub>e</sub> = 74, E<sub>s</sub> = 500 eV; needed for '
              'Nixon et al. (equal energies)'],
             ['Protocol, §2.3', 'Claims better recoil peak with 3CWZ', 'Not confirmed for ' + NH3 + ' at 74 eV; present it as '
              'a question the protocol tests'],
             ['Protocol, ref. [26] and §5', 'Moccia J. Chem. Phys. 40, 2164 (paper I) for both molecules; "Tables IV and IX"',
              'Cite 2164 for ' + CH4 + ' and 2176 (paper II, Table I) for ' + NH3 + ', as Q3 does'],
             ['Protocol, after Eq. (3)', 'r<sub>j</sub> called the "ejected electron"', 'r<sub>j</sub> are the target electrons'],
             ['Protocol, objective 2', '"l\'orbitale de NH₃" (singular)', 'The three valence orbitals are treated'],
             ['Protocol, §5', 'No numerical validation step', 'Add a grid-convergence step (section 4.3)']],
            [3.0 * cm, 5.6 * cm, 7.6 * cm])]

# ---- next steps
s += [P('6. Remaining work, ordered by cost and mapped to the protocol', h2),
      table([['#', 'Task', 'Protocol item', 'Cost'],
             ['1', 'PCI switch (C<sub>12</sub> = 1) and ' + NH3 + ' runs at 74 eV', 'O3', 'small'],
             ['2', 'Same runs with constant Z = 1 (3C/BBK limit) to isolate the effect of Z(r)', 'Expected result 1', 'small'],
             ['3', 'Digitize the 37 eV data of El Mir et al. and compare (runs already done)', 'O4', 'small'],
             ['4', CH4 + ' runs with the same code ('+ t2 + ', 2a<sub>1</sub>) vs Lahmam-Bennani et al. 2009 [35]', 'O2, O4', 'medium'],
             ['5', 'Distorted wave for the ejected electron in the Q3 potential', 'Expected result 2, perspectives', 'large'],
             ['6', 'Grid convergence / K-adapted grids, then integrated cross sections vs [38]', 'O4, method', 'medium'],
             ['7', 'Exchange amplitude, then Nixon et al. 2013 [39] (low energy)', 'O4', 'large']],
            [0.8 * cm, 9.6 * cm, 3.6 * cm, 2.2 * cm])]

s += [P('7. Files (repository dramoros95/dd, branch claude/read-and-explain-k311k1)', h2),
      table([['Path', 'Content'],
             ['nh3/nh3_opt_claude.for', 'NH<sub>3</sub> 3CWZ code (options isr, iorth)'],
             ['nh3/elmir_conditions/', 'Runs, experimental data, figures and PDF discussion for El Mir et al.'],
             ['nh3/theory_orthogonalization.md', 'Theory of the orthogonalization test'],
             ['nh3/theory_distorted_wave.md, nh3/diagnostics/', 'Ejected-wave diagnosis and proposed fix'],
             ['nh3/rezkallah_conditions/', 'EMS and DDCS tests, data, scripts, README'],
             ['nh3/progress/', 'This synthesis, Z(r) data and figure']], [6.5 * cm, 9.7 * cm])]

def footer(c, d):
    c.saveState(); c.setFont('DV', 7.5); c.setFillColor(colors.HexColor('#8a8985'))
    c.drawRightString(A4[0] - 2 * cm, 1.2 * cm, 'page %d' % d.page)
    c.drawString(2 * cm, 1.2 * cm, '3CWZ CH4/NH3: progress against the protocol')
    c.restoreState()

SimpleDocTemplate(OUT, pagesize=A4, leftMargin=2 * cm, rightMargin=2 * cm, topMargin=1.6 * cm, bottomMargin=1.8 * cm,
                  title='3CWZ CH4/NH3: progress against the protocol', author='').build(s, onFirstPage=footer, onLaterPages=footer)
print(OUT)
