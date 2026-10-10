"""The single introductory figure: physical query -> sources -> endpoint -> tests.

All labels, coordinates and diagram topology are native SVG. The six embedded
source icons, solvent halo and endpoint lens are conceptual illustration
components, not chemical structures, calculated fields or experimental data.
"""
import base64

from journal_vector_art import Art, ROOT, INK, MID, GRID, BLUE, TEAL, AMBER, COLORS, COUNTS
from journal_diagram import wire, heading

OUT = ROOT / 'paper/figures/journal'
ASSET = OUT / 'teaser_components'


def embed_component(c, name, x, y, width, height):
    """Embed one publication-ready raster; the selected component is pure white."""
    encoded = base64.b64encode((ASSET / f'{name}_white.png').read_bytes()).decode()
    c.add(f'<image x="{x}" y="{y}" width="{width}" height="{height}" '
          f'xlink:href="data:image/png;base64,{encoded}"/>')


def ribbon_texture(c, j, shape, source_y, model_y, color):
    """Clip one embedded illustration to one exact band, with no tiled seams."""
    encoded = base64.b64encode((ASSET / f'paper_texture_{j+1}.png').read_bytes()).decode()
    top = min(source_y, model_y) - 12
    height = abs(source_y - model_y) + 24
    clip_id = f'source-ribbon-clip-{j}'
    c.defs.append(f'<clipPath id="{clip_id}" clipPathUnits="userSpaceOnUse">'
                  f'<path d="{shape}"/></clipPath>')
    c.path(shape, color, 0, color)
    # The image is deliberately larger than the clip on all four sides.
    # An opaque, single image avoids antialiased seams from pattern repeats.
    c.add(f'<image x="249" y="{top-1}" width="192" height="{height+2}" '
          f'preserveAspectRatio="none" clip-path="url(#{clip_id})" '
          f'xlink:href="data:image/png;base64,{encoded}"/>')
    c.path(shape, color, .55)


def output(c, x, y, size=25):
    c.text(x, y, 'ΔG', size, TEAL, 500)
    c.text(x + size * 1.58, y + size * .22, 'hyd', size * .43, TEAL)


def structure_features(c, base_y=596):
    """Schematic structure featurization; the molecule is the verified NMA example."""
    offset = base_y - 596
    c.molecule(36, base_y, 43, 26, hydrogens=False)
    wire(c, (80, 609 + offset), (91, 609 + offset), MID, width=1)
    # Bit and property glyphs describe feature types, not actual NMA values.
    bits = ('110010101101', '011001011010', '101100101001')
    for row, pattern in enumerate(bits):
        for col, bit in enumerate(pattern):
            c.rect(96 + col * 4.9, 599 + offset + row * 6.2, 3.8, 4.7,
                   INK if bit == '1' else '#E3EBEE', radius=.45)
    c.text(162, 613 + offset, '+', 15, MID, 400, 'middle')
    for j, height in enumerate((8, 15, 11, 18, 10, 16, 13)):
        c.rect(175 + j * 7.7, 618 + offset - height, 5.2, height,
               MID if j % 2 else '#9AABB3', radius=.6)
    c.text(57, 628 + offset, 'SMILES', 9.2, MID, 500, 'middle')
    c.text(124, 628 + offset, 'fingerprint', 9.2, MID, 500, 'middle')
    c.text(201, 628 + offset, 'RDKit features', 9.2, MID, 500, 'middle')


def molecule_records(c, x, y, columns, rows, first_count, first_color,
                     second_color=None, step_x=7.1, step_y=5.2):
    """Exact-count cohort mosaic; tiles stand for records, not molecular identity."""
    for k in range(columns * rows):
        col = k % columns
        row = k // columns
        color = first_color if k < first_count else second_color
        assert color is not None
        c.rect(x + col * step_x, y + row * step_y,
               step_x - 1.6, step_y - 1.4, color, radius=.6)


def overview(data, preview):
    """One non-repeating overview of source learning, student fit and evaluation."""
    c = Art(860, 'SolvAI: reusable solvent information from distinct sources')

    # a | Show the motivation exactly once. N-methylacetamide connectivity is
    # verified in Art; the water shell is explicitly labelled schematic.
    heading(c, 10, 28, 'a', 'The cost of a new solvent calculation')
    c.molecule(28, 61, 80, 61)
    wire(c, (112, 94), (145, 94), BLUE, width=1.7)
    embed_component(c, 'solvent_halo', 151, 36, 118, 118)
    # The solute is a native, connectivity-checked NMA rendering. The blue halo
    # is a conceptual illustration, not a simulated water configuration.
    c.molecule(181, 78, 56, 38)
    wire(c, (272, 94), (300, 94), BLUE, width=1.7)
    output(c, 313, 101, 19)
    c.text(33, 163, 'New molecule', 11, MID, 400)
    c.text(162, 163, 'Sample water · schematic', 11, MID, 400)
    c.text(424, 79, 'A new physical calculation is costly', 13, INK, 600)
    c.text(424, 100, 'when it must be repeated for every query.', 11.4, MID)
    c.text(424, 122, 'SolvAI learns from existing solvent information.', 11.4, TEAL, 500)
    c.line(10, 184, 710, 184, GRID, 1)

    # b | Each distinct source dataset teaches its own structure-to-response
    # model group. Counts and maps are exact, while the icons are conceptual.
    heading(c, 10, 207, 'b', 'Six source datasets teach reusable response maps')
    c.text(35, 227, 'Different targets; one common molecular input', 11.1, MID)
    c.text(80, 249, 'SOURCE DATA AND TARGET', 9.4, MID, 600)
    c.text(451, 249, 'SOURCE MODEL GROUP', 9.4, MID, 600, 'middle')
    c.text(574, 249, 'PREDICTED RESPONSE', 9.4, MID, 600, 'middle')
    names = ['CombiSolv-QM', 'SoluteML', 'OpenFF', 'GBn2', 'MolSolv', 'ConfSolv']
    targets = ['COSMOtherm water', 'Abraham interaction axes',
               'Explicit-water + correction', 'Implicit-water + correction',
               'SMD(water)', 'Conformer-response summaries']
    for j, (name, target, color, n) in enumerate(zip(names, targets, COLORS, COUNTS)):
        y = 260 + j * 44
        icon = base64.b64encode((ASSET / f'source_icon_{j+1}.png').read_bytes()).decode()
        c.add(f'<image x="25" y="{y-3}" width="44" height="44" '
              f'xlink:href="data:image/png;base64,{icon}"/>')
        c.text(80, y + 16, name, 11.8, color, 600)
        c.text(80, y + 32, target, 9.5, MID)
        source_y = y + 21
        # The six streams contract toward their common middle as they approach
        # the learned maps. Each band also tapers smoothly: no oscillatory wave.
        model_y = 391 + .82 * (source_y - 391)
        # Generated color lives only inside a deterministic vector silhouette.
        # The taper carries no scientific magnitude or measured field.
        shape = (f'M250,{source_y-10} C318,{source_y-10} 372,{model_y-5.5} 439,{model_y-5.5} '
                 f'L439,{model_y+5.5} C372,{model_y+5.5} 318,{source_y+10} '
                 f'250,{source_y+10} Z')
        ribbon_texture(c, j, shape, source_y, model_y, color)
        c.circle(453, model_y, 14, 'white', color, 1.15)
        c.text(453, model_y + 4, 'φ' + str(j+1), 11, color, 500, 'middle')
        wire(c, (467, model_y), (482, model_y), color, width=1.15)
        for k in range(n):
            c.rect(489 + k * 14, model_y - 8, 10, 16, color, radius=.8)
        c.text(608, model_y + 4, f'{n} coordinate' + ('' if n == 1 else 's'), 10, MID)
    c.text(359, 539, 'Freeze six model groups  →  15 predicted response coordinates',
           11.4, TEAL, 600, 'middle')
    c.line(10, 552, 710, 552, GRID, 1)

    # c | Continue the same 15 colored coordinates from b, now alongside the
    # independent SMILES-derived structure representation. The central lens is
    # an architecture-neutral learned endpoint, not a depicted network layer.
    heading(c, 10, 576, 'c', 'Learn hydration from experimental labels')
    c.text(35, 612, 'r̂(x)  ·  predicted solvent responses', 10.3, TEAL, 600)
    c.responses(35, 621, width=201, height=17)
    c.text(241, 635, '15', 10.3, TEAL, 600)
    c.text(35, 657, 's(x)  ·  structure from SMILES', 10.3, MID, 600)
    structure_features(c, base_y=660)
    # The two streams contract into a shared learned mapping.
    c.path('M258,630 C295,630 310,660 352,660', TEAL, 1.7)
    c.path('M258,678 C295,678 312,671 352,671', MID, 1.7)
    c.path('M352,660 Q363,665 352,671', TEAL, 1.2)
    # The generated glass motif intentionally has no tree, graph or neural
    # topology: all tested endpoint families can occupy this position.
    embed_component(c, 'endpoint_lens', 374, 625, 132, 88)
    wire(c, (354, 665), (393, 665), TEAL, width=1.6)
    c.text(440, 675, 'h', 27, TEAL, 500, 'middle', italic=True)
    c.text(440, 707, 'Endpoint learner', 10.8, TEAL, 600, 'middle')
    c.text(440, 612, 'Experimental ΔG labels', 10.3, AMBER, 600, 'middle')
    wire(c, (440, 616), (440, 643), AMBER, width=1.1)
    wire(c, (486, 665), (600, 665), TEAL, width=1.5)
    output(c, 615, 673, 23)
    c.text(35, 726, 'Inference: SMILES only; no new solvent simulation or measurement.',
           10.8, TEAL, 500)
    c.line(10, 738, 710, 738, GRID, 1)

    # d | The two exclusions answer different questions. The strict 97 are a
    # subset of 220, so the right mosaic is a partition, not a third cohort.
    heading(c, 10, 764, 'd', 'Evaluation separates endpoint and source exposure')
    c.line(352, 781, 352, 856, GRID, 1)
    c.text(30, 798, 'ARROW-85', 13.2, BLUE, 600)
    molecule_records(c, 31, 808, 17, 5, 85, BLUE, step_x=7.6, step_y=8.1)
    c.text(174, 817, '85 absent from all six', 10.4, INK, 500)
    c.text(174, 833, 'supervised source sets', 10.4, INK)
    c.text(174, 850, 'Endpoint: out of fold', 10.4, MID)
    c.text(369, 798, 'External-220', 13.2, BLUE, 600)
    molecule_records(c, 370, 806, 22, 10, 123, MID, TEAL,
                     step_x=7.0, step_y=5.2)
    c.text(540, 818, '220 endpoint-disjoint', 10.1, INK, 500)
    c.rect(540, 827, 7, 7, MID, radius=.5)
    c.text(551, 834, '123 source-exposed', 9.7, MID, 500)
    c.rect(540, 843, 7, 7, TEAL, radius=.5)
    c.text(551, 850, '97 source-disjoint', 9.7, TEAL, 600)
    return c.save(OUT / 'F1_overview', preview / 'F1_overview.png' if preview else None)
