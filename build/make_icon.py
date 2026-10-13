"""Draws the EasyCGTree icon (a cladogram inside the two strands of a genome with core genes, on a navy
tile) and writes the icon files:
icon.png (1024 px), icon.ico (Windows), icon.icns (macOS), and icon_b64.txt (256 px PNG as base64, which
is pasted into EasyCGTree_GUI.py as the window icon). Needs Pillow:  python3 make_icon.py"""
import base64
import io
from pathlib import Path

from PIL import Image, ImageDraw

S = 4                      # supersampling
N = 1024 * S
NAVY_TOP, NAVY_BOTTOM = (36, 92, 140), (19, 52, 84)
LINE = (255, 255, 255)
TIP = (110, 196, 184)      # muted teal
AMBER = (232, 176, 92)


def draw():
    """Two rings (the two strands of a genome) with white blocks (core genes) around a cladogram."""
    grad = Image.new("RGB", (1, N))
    for y in range(N):
        t = y / (N - 1)
        grad.putpixel((0, y), tuple(round(a + (b - a) * t) for a, b in zip(NAVY_TOP, NAVY_BOTTOM)))
    grad = grad.resize((N, N))
    mask = Image.new("L", (N, N), 0)
    m = 100 * S
    ImageDraw.Draw(mask).rounded_rectangle([m, m, N - m, N - m], radius=185 * S, fill=255)
    img = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    img.paste(grad, (0, 0), mask)
    d = ImageDraw.Draw(img)

    def p(v):
        return v * S
    cx, cy = 512, 512
    R1, R2, W = 330, 282, 26                 # outer and inner strand
    for R, col in ((R1, TIP), (R2, AMBER)):
        d.ellipse([p(cx - R), p(cy - R), p(cx + R), p(cy + R)], outline=col, width=p(W))
    genes = {R1: ((10, 30), (70, 88), (150, 172), (232, 246), (290, 312)),
             R2: ((40, 58), (110, 128), (196, 214), (262, 278), (332, 350))}
    for R, arcs in genes.items():
        box = [p(cx - R - 7), p(cy - R - 7), p(cx + R + 7), p(cy + R + 7)]
        for a0, a1 in arcs:
            d.arc(box, a0, a1, fill=LINE, width=p(W + 14))
    w = p(30)

    def seg(x1, y1, x2, y2):
        d.line([(p(x1), p(y1)), (p(x2), p(y2))], fill=LINE, width=w)
        for x, y in ((x1, y1), (x2, y2)):
            d.ellipse([p(x) - w / 2, p(y) - w / 2, p(x) + w / 2, p(y) + w / 2], fill=LINE)
    # topology ((A,B),(C,(D,E))), root on the left
    tips = [372, 442, 512, 582, 652]
    tipx = 660
    n1 = (530, (tips[0] + tips[1]) / 2)
    n3 = (560, (tips[3] + tips[4]) / 2)
    n2 = (450, (tips[2] + n3[1]) / 2)
    root = (355, (n1[1] + n2[1]) / 2)
    seg(300, root[1], root[0], root[1])
    seg(root[0], n1[1], root[0], n2[1])
    seg(root[0], n1[1], n1[0], n1[1])
    seg(root[0], n2[1], n2[0], n2[1])
    seg(n1[0], tips[0], n1[0], tips[1])
    seg(n2[0], tips[2], n2[0], n3[1])
    seg(n2[0], n3[1], n3[0], n3[1])
    seg(n3[0], tips[3], n3[0], tips[4])
    for y, x0 in zip(tips, (n1[0], n1[0], n2[0], n3[0], n3[0])):
        seg(x0, y, tipx, y)
    r = p(22)
    for y in tips:
        d.ellipse([p(tipx + 16) - r, p(y) - r, p(tipx + 16) + r, p(y) + r], fill=TIP, outline=LINE, width=p(7))
    return img.resize((1024, 1024), Image.LANCZOS)


def main():
    here = Path(__file__).resolve().parent
    img = draw()
    img.save(here / "icon.png")
    img.save(here / "icon.ico", sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
    img.save(here / "icon.icns")
    buf = io.BytesIO()
    img.resize((256, 256), Image.LANCZOS).save(buf, "PNG", optimize=True)
    (here / "icon_b64.txt").write_text(base64.b64encode(buf.getvalue()).decode("ascii"))
    print("written: icon.png, icon.ico, icon.icns, icon_b64.txt (%d bytes)" % len(buf.getvalue()))


if __name__ == "__main__":
    main()
