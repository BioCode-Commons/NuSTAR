"""Export the panel walkthrough from browser-captured, intrinsic-size PNGs.

Capture implementation-comparison.html?figure=flow&capture=1 and
?figure=diff&capture=1 with a full-page browser screenshot first. Store them as
docs/assets/implementation-flow.png and implementation-differences.png.
Then run this script from any directory. Pillow is needed only for this export.
"""
from pathlib import Path
import xml.etree.ElementTree as ET
from PIL import Image, ImageDraw

assets = Path(__file__).resolve().parents[2] / 'docs/assets'
namespace = {'svg': 'http://www.w3.org/2000/svg'}
for name in ('implementation-flow', 'implementation-differences'):
    svg = ET.parse(assets / (name + '.svg')).getroot()
    width, height = int(svg.attrib['width']), int(svg.attrib['height'])
    image = Image.open(assets / (name + '.png')).convert('RGB')
    if image.width < width or image.height != height:
        raise ValueError(f'{name}: screenshot must show the full SVG at intrinsic size')
    image = image.crop((0, 0, width, height))
    image.save(assets / (name + '.png'))
    if name == 'implementation-differences':
        frames = []
        for group in svg.findall('svg:g', namespace):
            if 'data-difference' not in group.attrib:
                continue
            panel = group.find('svg:rect', namespace)
            x, y, w, h = (int(panel.attrib[key]) for key in ('x', 'y', 'width', 'height'))
            frame = image.copy()
            draw = ImageDraw.Draw(frame)
            draw.rounded_rectangle((x, y, x + w, y + h), radius=8, outline='#58a6ff', width=3)
            draw.rectangle((x, y + 8, x + 4, y + h - 8), fill='#58a6ff')
            frames.append(frame.quantize(colors=128))
        if not frames:
            raise ValueError('No difference panels found')
        frames[0].save(assets / 'implementation-differences.gif', save_all=True,
                       append_images=frames[1:], duration=1800, loop=0, optimize=True)
        print(f'Exported {len(frames)} animated difference panels')
