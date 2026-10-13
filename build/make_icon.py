"""Draws the EasyCGTree icon (a cladogram on a light disc inside the two strands of a genome with core
genes) and writes the icon files:
icon.png (1024 px), icon.ico (Windows), icon.icns (macOS), and icon_b64.txt (256 px PNG as base64, which
is pasted into EasyCGTree_GUI.py as the window icon). Needs Pillow:  python3 make_icon.py"""
import base64
import io
from pathlib import Path

from PIL import Image, ImageDraw

S = 4                      # supersampling
N = 1024 * S
NAVY = (31, 78, 121)
DISC = (244, 247, 250)     # light disc
BLOCK = 1 / 3              # length of the core-gene blocks (fraction of the first design)
LINE = (255, 255, 255)
TIP = (110, 196, 184)      # muted teal
AMBER = (232, 176, 92)


def draw():
    """A light disc inside two rings (the two strands of a genome) with short navy blocks (core genes),
    and a navy cladogram on the disc. No background tile: it works on light and dark backgrounds."""
    img = Image.new("RGBA", (N, N), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)

    def p(v):
        return v * S
    cx = cy = 512
    R1, R2, W = 420, 362, 34                 # outer and inner strand
    R0 = R2 - W / 2 - 6
    d.ellipse([p(cx - R0), p(cy - R0), p(cx + R0), p(cy + R0)], fill=DISC)
    for R, col in ((R1, TIP), (R2, AMBER)):
        d.ellipse([p(cx - R), p(cy - R), p(cx + R), p(cy + R)], outline=col, width=p(W))
    genes = {R1: ((10, 30), (70, 88), (150, 172), (232, 246), (290, 312)),
             R2: ((40, 58), (110, 128), (196, 214), (262, 278), (332, 350))}
    for R, arcs in genes.items():
        box = [p(cx - R - 8), p(cy - R - 8), p(cx + R + 8), p(cy + R + 8)]
        for a0, a1 in arcs:
            m, h = (a0 + a1) / 2, (a1 - a0) * BLOCK / 2
            d.arc(box, m - h, m + h, fill=NAVY, width=p(W + 16))
    w = p(38)
    k = 1.28                                   # the tree is drawn on a 1024 grid and scaled around the centre

    def seg(x1, y1, x2, y2):
        x1, y1, x2, y2 = cx + (x1 - 512) * k, cy + (y1 - 512) * k, cx + (x2 - 512) * k, cy + (y2 - 512) * k
        d.line([(p(x1), p(y1)), (p(x2), p(y2))], fill=NAVY, width=w)
        for x, y in ((x1, y1), (x2, y2)):
            d.ellipse([p(x) - w / 2, p(y) - w / 2, p(x) + w / 2, p(y) + w / 2], fill=NAVY)
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
    r = p(28)
    for y in tips:
        X, Y = cx + (tipx + 16 - 512) * k, cy + (y - 512) * k
        d.ellipse([p(X) - r, p(Y) - r, p(X) + r, p(Y) + r], fill=TIP, outline=LINE, width=p(8))
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
