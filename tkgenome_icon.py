#!/usr/bin/env python3
"""Icon and macOS app-bundle builder for tkgenome (standard library only).

  tkgenome_icon.py png SIZE [OUT.png]   write the icon as a PNG
  tkgenome_icon.py b64                  print the 128px PNG as base64 (the
                                        ICON_PNG constant embedded in tkgenome)
  tkgenome_icon.py app [DEST_DIR]       build tkgenome.app (macOS) next to
                                        this script, or in DEST_DIR

The icon is a miniature of the map: a dark tile with a circular genome,
forward-strand arcs outside, reverse-strand arcs inside, and one highlighted
feature in red.
"""

import base64
import math
import os
import plistlib
import shutil
import struct
import sys
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))

BG = (27, 42, 65)
BACKBONE = (150, 160, 175)
PLUS = (47, 111, 208)
MINUS = (232, 137, 43)
GREEN = (46, 158, 72)
RED = (255, 31, 31)
WHITE = (240, 240, 240)

PLUS_ARCS = [(8, 30), (34, 38), (66, 100), (106, 140), (150, 200),
             (206, 214), (250, 300), (310, 350)]
MINUS_ARCS = [(0, 22), (28, 70), (80, 120), (128, 134), (160, 220),
              (228, 262), (270, 330), (338, 356)]
TRNA_ARCS = [(100, 104), (224, 227)]
HIGHLIGHT = (40, 62)


def _in(theta, arcs):
    return any(a <= theta < b for a, b in arcs)


def sample(x, y):
    """RGBA of the icon at (x, y) in the unit square, or None (transparent)."""
    # rounded-square tile
    cx = min(max(x, 0.2), 0.8)
    cy = min(max(y, 0.2), 0.8)
    if math.hypot(x - cx, y - cy) > 0.2 or not (0 <= x <= 1 and 0 <= y <= 1):
        return None
    dx, dy = x - 0.5, y - 0.5
    r = math.hypot(dx, dy)
    th = math.degrees(math.atan2(dx, -dy)) % 360.0
    if 0.325 <= r <= 0.415 and _in(th, [HIGHLIGHT]):
        return RED
    if 0.34 <= r <= 0.40 and _in(th, PLUS_ARCS):
        return PLUS
    if 0.20 <= r <= 0.26 and _in(th, MINUS_ARCS):
        return MINUS
    if 0.34 <= r <= 0.40 and _in(th, TRNA_ARCS):
        return GREEN
    if 0.20 <= r <= 0.26 and _in(th, [(100, 104)]):
        return GREEN
    if abs(r - 0.30) <= 0.006:
        return BACKBONE
    if 0.425 <= r <= 0.46 and (th < 1.5 or th > 358.5):
        return WHITE
    return BG


def render_png(size):
    n = 4 if size <= 256 else 2
    rows = []
    for py in range(size):
        row = bytearray([0])                     # PNG filter type 0
        for px in range(size):
            ra = ga = ba = aa = 0
            for sy in range(n):
                for sx in range(n):
                    c = sample((px + (sx + 0.5) / n) / size,
                               (py + (sy + 0.5) / n) / size)
                    if c:
                        ra += c[0]
                        ga += c[1]
                        ba += c[2]
                        aa += 1
            if aa:
                row += bytes((ra // aa, ga // aa, ba // aa,
                              255 * aa // (n * n)))
            else:
                row += b"\0\0\0\0"
        rows.append(bytes(row))

    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data +
                struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff))
    return (b"\x89PNG\r\n\x1a\n" +
            chunk(b"IHDR", struct.pack(">IIBBBBB", size, size, 8, 6, 0, 0, 0)) +
            chunk(b"IDAT", zlib.compress(b"".join(rows), 9)) +
            chunk(b"IEND", b""))


# ICNS entry types holding PNG data, with the pixel size of each
ICNS_TYPES = [(b"icp4", 16), (b"icp5", 32), (b"ic11", 32), (b"ic12", 64),
              (b"ic07", 128), (b"ic13", 256), (b"ic08", 256), (b"ic14", 512),
              (b"ic09", 512)]


def make_icns():
    """An .icns file written directly (no iconutil needed)."""
    pngs = {}
    body = b""
    for tag, size in ICNS_TYPES:
        if size not in pngs:
            pngs[size] = render_png(size)
        body += tag + struct.pack(">I", len(pngs[size]) + 8) + pngs[size]
    return b"icns" + struct.pack(">I", len(body) + 8) + body


def build_app(dest):
    app = os.path.join(dest, "tkgenome.app")
    if os.path.exists(app):
        shutil.rmtree(app)
    macos = os.path.join(app, "Contents", "MacOS")
    res = os.path.join(app, "Contents", "Resources")
    os.makedirs(macos)
    os.makedirs(res)
    shutil.copy(os.path.join(HERE, "tkgenome"),
                os.path.join(res, "tkgenome.py"))

    launcher = os.path.join(macos, "tkgenome")
    with open(launcher, "w") as fh:
        fh.write('#!/bin/sh\nexec "%s" "$(dirname "$0")/../Resources/'
                 'tkgenome.py" "$@"\n' % sys.executable)
    os.chmod(launcher, 0o755)

    with open(os.path.join(res, "tkgenome.icns"), "wb") as fh:
        fh.write(make_icns())

    exts = ["gb", "gbk", "gbff", "genbank"]
    with open(os.path.join(app, "Contents", "Info.plist"), "wb") as fh:
        plistlib.dump({
            "CFBundleName": "tkgenome",
            "CFBundleDisplayName": "tkgenome",
            "CFBundleIdentifier": "org.badgerjh.tkgenome",
            "CFBundleExecutable": "tkgenome",
            "CFBundleIconFile": "tkgenome",
            "CFBundlePackageType": "APPL",
            "CFBundleShortVersionString": "1.0",
            "NSHighResolutionCapable": True,
            "CFBundleDocumentTypes": [{
                "CFBundleTypeName": "GenBank file",
                "CFBundleTypeRole": "Viewer",
                "CFBundleTypeExtensions": exts,
            }],
        }, fh)
    return app


def main():
    args = sys.argv[1:]
    if args[:1] == ["png"] and len(args) >= 2:
        data = render_png(int(args[1]))
        with open(args[2] if len(args) > 2 else "tkgenome.png", "wb") as fh:
            fh.write(data)
    elif args[:1] == ["b64"]:
        b = base64.b64encode(render_png(128)).decode()
        print("\n".join(b[i:i + 76] for i in range(0, len(b), 76)))
    elif args[:1] == ["app"]:
        print(build_app(args[1] if len(args) > 1 else HERE))
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
