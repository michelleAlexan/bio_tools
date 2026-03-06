#%%
import re

input_file = "/Users/michellealexander/Downloads/tree-31.svg"
output_file = "/Users/michellealexander/Downloads/tree-31_fix.svg"



with open(input_file, "r", encoding="utf-8") as f:
    svg = f.read()
    print(svg)

# # <path> elements 
# svg = re.sub(r'<path[^>]*>', '', svg)
def make_path_transparent(match):
    tag = match.group(0)

    # Remove existing fill and stroke definitions
    tag = re.sub(r'fill="[^"]*"', '', tag)
    tag = re.sub(r'stroke="[^"]*"', '', tag)

    # Add transparent styling
    tag = tag.replace('<path', '<path fill="darkgrey" stroke="darkgrey"')

    return tag

svg = re.sub(r'<path[^>]*>', make_path_transparent, svg)


# 2️⃣ Recolor all <line> elements (branches) to darkgrey
def recolor_line(match):
    line_tag = match.group(0)
    # Replace stroke if exists, else add it
    if 'stroke=' in line_tag:
        line_tag = re.sub(r'stroke="[^"]*"', 'stroke="darkgrey"', line_tag)
    else:
        line_tag = line_tag.replace('<line', '<line stroke="darkgrey"')
    # Optional: ensure a default stroke-width
    if 'stroke-width=' not in line_tag:
        line_tag = line_tag.replace('<line', '<line stroke-width="1"')
    return line_tag

svg = re.sub(r'<line[^>]*>', recolor_line, svg)


# 2️⃣ Find all rects
# 1️⃣ Remove all nodebox rectangles WITHOUT fill
# 1️⃣ Remove all rects that have class="nodebox"
svg = re.sub(r'<rect[^>]*class="nodebox"[^>]*/?>', '', svg)

# 2️⃣ Make all <path> white (so they don't cover things)
def make_path_white(match):
    tag = match.group(0)
    tag = re.sub(r'fill="[^"]*"', '', tag)
    tag = re.sub(r'stroke="[^"]*"', '', tag)
    return tag.replace('<path', '<path fill="white" stroke="none"')

# svg = re.sub(r'<path[^>]*>', make_path_white, svg)

# 3️⃣ Recolor all <line> elements to darkgrey
def recolor_line(match):
    line_tag = match.group(0)

    if 'stroke=' in line_tag:
        line_tag = re.sub(r'stroke="[^"]*"', 'stroke="darkgrey"', line_tag)
    else:
        line_tag = line_tag.replace('<line', '<line stroke="darkgrey"')

    if 'stroke-width=' not in line_tag:
        line_tag = line_tag.replace('<line', '<line stroke-width="1"')

    return line_tag

svg = re.sub(r'<line[^>]*>', recolor_line, svg)


svg = re.sub(
    r'<path\b[^>]*class="collapsed"[^>]*>',
    lambda m: (
        re.sub(r'fill="[^"]*"', '', m.group(0))
        .replace('<path', '<path fill="none" stroke="none"')
    ),
    svg
)

# 4️⃣ Move remaining rects to top layer
rect_pattern = r'<rect[^>]*>'
remaining_rects = re.findall(rect_pattern, svg)

# remove them from current positions
svg = re.sub(rect_pattern, '', svg)

# append at end (top layer)
svg = svg.replace('</svg>', ''.join(remaining_rects) + '\n</svg>')

with open(output_file, "w", encoding="utf-8") as f:
    f.write(svg)
print("SVG fixed! Branch lines are darkgrey, nodeboxes intact, paths removed.")

# %%


with open(input_file, "r", encoding="utf-8") as f:
    svg = f.read()


# 1️⃣ Remove only nodeboxes WITHOUT fill
def remove_uncolored_nodeboxes(match):
    tag = match.group(0)

    # Check for fill attribute (fill="something")
    has_fill_attr = re.search(r'fill\s*=\s*"[^"]+"', tag, flags=re.IGNORECASE)

    # Check for fill inside style attribute (style="...fill: ...;...")
    has_fill_style = re.search(r'style\s*=\s*"[^"]*fill\s*:\s*[^;"]+', tag, flags=re.IGNORECASE)

    if has_fill_attr or has_fill_style:
        return tag  # keep colored boxes

    return ''  # remove background ones

svg = re.sub(r'<rect[^>]*class="nodebox"[^>]*/?>', remove_uncolored_nodeboxes, svg, flags=re.IGNORECASE)

# 2️⃣ Make all <path> invisible (fill white, stroke none)
def make_path_white(match):
    tag = match.group(0)
    tag = re.sub(r'fill="[^"]*"', '', tag)
    tag = re.sub(r'stroke="[^"]*"', '', tag)
    return tag.replace('<path', '<path fill="white" stroke="none"')

# svg = re.sub(r'<path[^>]*>', make_path_white, svg, flags=re.IGNORECASE)

# 3️⃣ Recolor <line> elements to darkgrey
def recolor_line(match):
    tag = match.group(0)

    if 'stroke=' in tag:
        tag = re.sub(r'stroke="[^"]*"', 'stroke="darkgrey"', tag)
    else:
        tag = tag.replace('<line', '<line stroke="darkgrey"')

    if 'stroke-width=' not in tag:
        tag = tag.replace('<line', '<line stroke-width="1"')

    return tag

svg = re.sub(r'<line[^>]*>', recolor_line, svg, flags=re.IGNORECASE)

# 4️⃣ Optional: move colored rects to top layer
rect_pattern = r'<rect[^>]*class="nodebox"[^>]*/?>'
colored_rects = re.findall(rect_pattern, svg, flags=re.IGNORECASE)
svg = re.sub(rect_pattern, '', svg, flags=re.IGNORECASE)
svg = svg.replace('</svg>', ''.join(colored_rects) + '\n</svg>')

with open(output_file, "w", encoding="utf-8") as f:
    f.write(svg)

print("✅ SVG processed! Only colored nodeboxes remain, paths are white, lines are darkgrey.")


# %%
import re

input_file = "/Users/michellealexander/Downloads/tree_fix1.svg"

with open(input_file, "r", encoding="utf-8") as f:
    svg = f.read()

# Find all colored nodeboxes
rect_pattern = r'<rect[^>]*class="nodebox"[^>]*/?>'
colored_rects = re.findall(rect_pattern, svg, flags=re.IGNORECASE)

legend_data = []

for rect in colored_rects:
    # Extract color from style
    style_match = re.search(r'style="[^"]*fill\s*:\s*([^;"]+)', rect)
    color = style_match.group(1) if style_match else 'none'
    
    # Extract the data-info attribute for label
    info_match = re.search(r'data-info="([^"]+)"', rect)
    info = info_match.group(1) if info_match else ''
    
    # Optional: clean HTML entities
    info = info.replace('&lt;', '<').replace('&gt;', '>').replace('&amp;', '&')
    
    legend_data.append((color, info.split(" ")[-1]))

# Print legend info
for color, info in legend_data:
    print(f"{color} -> {info}")


output_file_legend = "/Users/michellealexander/Downloads/tree_legend.svg"


# 2️⃣ Generate a separate SVG for the legend
legend_width = 300
legend_height = 50 + 30 * len(legend_data)
box_size = 20
spacing = 30

legend_svg = f'<svg xmlns="http://www.w3.org/2000/svg" width="{legend_width}" height="{legend_height}">\n'

for i, (color, info) in enumerate(legend_data):
    y_pos = 20 + i * spacing
    # Colored rectangle with opacity 0.4
    legend_svg += f'  <rect x="10" y="{y_pos}" width="{box_size}" height="{box_size}" fill="{color}" fill-opacity="0.4" stroke="black"/>\n'
    # Text label
    legend_svg += f'  <text x="{10 + box_size + 5}" y="{y_pos + box_size*0.75}" font-size="16">{info}</text>\n'

legend_svg += '</svg>'

with open(output_file_legend, "w", encoding="utf-8") as f:
    f.write(legend_svg)
# %%



# %%----------------##THIS WORKS FOR CIRCULAR TREEEEES!!------------------------------

import re

input_file = "/Users/michellealexander/Downloads/tree-36.svg"
output_file = "/Users/michellealexander/Downloads/tree-36_fix.svg"

with open(input_file, "r", encoding="utf-8") as f:
    svg = f.read()

# --------------------------------------------------
# Handle <path> correctly
# --------------------------------------------------
def process_path(match):
    tag = match.group(0)

    # Extract class attribute
    class_match = re.search(r'class="([^"]*)"', tag)
    classes = class_match.group(1) if class_match else ""

    # 1️⃣ If this is a branch path → keep + recolor
    if "line" in classes:
        tag = re.sub(r'stroke="[^"]*"', '', tag)
        tag = re.sub(r'fill="[^"]*"', '', tag)
        return tag.replace('<path', '<path stroke="darkgrey" fill="none" stroke-width="1"')

    # 2️⃣ If collapsed overlay → hide
    if "collapsed" in classes:
        return tag.replace('<path', '<path fill="none" stroke="none" opacity="0"')

    # 3️⃣ Everything else → hide
    return tag.replace('<path', '<path fill="none" stroke="none" opacity="0"')


svg = re.sub(r'<path[^>]*>', process_path, svg)


# --------------------------------------------------
# 2️⃣ Recolor ONLY <line> elements (true branches)
# --------------------------------------------------
def recolor_line(match):
    tag = match.group(0)

    # replace or add stroke
    if 'stroke=' in tag:
        tag = re.sub(r'stroke="[^"]*"', 'stroke="darkgrey"', tag)
    else:
        tag = tag.replace('<line', '<line stroke="darkgrey"')

    # ensure stroke width
    if 'stroke-width=' not in tag:
        tag = tag.replace('<line', '<line stroke-width="1"')

    return tag

svg = re.sub(r'<line[^>]*>', recolor_line, svg)


# --------------------------------------------------
# 3️⃣ Remove nodebox rectangles (optional)
# --------------------------------------------------
svg = re.sub(r'<rect[^>]*class="nodebox"[^>]*/?>', '', svg)


# --------------------------------------------------
# 4️⃣ Move remaining rects to top layer
# --------------------------------------------------
rect_pattern = r'<rect[^>]*>'
remaining_rects = re.findall(rect_pattern, svg)

svg = re.sub(rect_pattern, '', svg)
svg = svg.replace('</svg>', ''.join(remaining_rects) + '\n</svg>')


# --------------------------------------------------
# Save
# --------------------------------------------------
with open(output_file, "w", encoding="utf-8") as f:
    f.write(svg)

print("SVG fixed correctly.")







# %%
import svgwrite

# -------------------------------
# 1️⃣ COLORS
# -------------------------------
GROUP_COLORS = {
    "Algae": "#574104",
    "Lycophytes": "#ab730c",
    "Liverworts": "#faaf00",
    "Ferns": "#c1d717",
    "Mosses": "#89be86",
    "Gymnosperms": "#076247",
    "Early Angiosperms": "#3BA0BC",
    "Monocots": "#7502d9",
    "Dicots": "#edc5ec",
}

COLORS_2ODD_FUNCTION ={
    "AOP2" : "#498BA0",
    "AOP3" : "#3D6F80",
    "DPS" : "#377672",
    "GA20ox" : "#6aa84f",
    "C20_GA2ox": "#406732",
    "C19_GA2ox": "#b6d7a8",
    "DAO" : "#539185",
    "GA3ox" : "#64b541",
    "GA13ox" : "#a0bd94",
    "GA7ox" : "#e8eed0",
    "2ODDC23" : "#fff2cc",
    "LFS" : "#f4e8c3",
    "2OG1" : "#ffe599",
    "C2H" : "#ffd966",
    "F6H" : "#ec8612",
    "S8H" : "#f57829",
    "GSLOH" : "#ec640f",
    "GRS" : "#e5ac00",
    "TIIAS" : "#af8300",
    "D4H" : "#bf7979",
    "BX6" : "#e0bbbb",
    "FNSI" : "#c492cc",
    "F3H" : "#f4cccc",
    "H6H" : "#e9d0db",
    "IDS" : "#cd87a6",
    "SLC" : "#cfe2f3",
    "GIM" : "#d0d2e5",
    "M2H" : "#c27ba0",
    "M2H_weak" : "#e688b8",
    "DMR6" : "#c2d3e2",
    "S5H" : "#abbbc9",
    "S3H" : "#95a3af",
    "FLS" : "#b4a7d6",
    "DAH": "#784fe1",
    "ANS" : "#8e7cc3",
    "JOX" : "#6fa8dc",
    "ACCO" : "#3d85c6",
    "T6OD" : "#3371a8",
    "COD" : "#316ca2",
    "SRG" : "#316a9f",
    "LBO" : "#2c6190",
    "NCS" : "#1a3852",
    "T2OGD": "#3B8274",

}

# -------------------------------
# 2️⃣ FUNCTION ORDER
# -------------------------------
leaf_function_order = [
    'GRS','GSLOH','D4H','TIIAS','BX6','S8H','F6H','C2H','GA7ox','2ODDC23','2OG1','LFS',
    'C19_GA2ox','GA3ox','GA13ox','GA20ox','C20_GA2ox','DAO','T2OGD', 'DPS','AOP2','AOP3',
    'JOX','LBO','NCS','M2H_weak','SRG','COD','T6OD','ACCO',
    'ANS','FLS','DAH','FNSI','F3H','S3H','S5H','DMR6','H6H','IDS','M2H','GIM','SLC'
]

# remove duplicates preserving order
seen = set()
ordered_functions = []
for f in leaf_function_order:
    if f not in seen:
        ordered_functions.append(f)
        seen.add(f)

# -------------------------------
# 3️⃣ CREATE SVG
# -------------------------------
dwg = svgwrite.Drawing("legend_lines.svg", size=("600px", "2000px"))

y = 30
x_box = 30
x_text = 80
box_size = 18
line_spacing = 28
branch_length = 40  # horizontal line length for functions

# ---- SECTION 1: Plant groups (squares) ----
dwg.add(dwg.text("Plant groups", insert=(30, y), font_size="20px", font_weight="bold"))
y += 40

for group, color in GROUP_COLORS.items():
    dwg.add(dwg.rect(
        insert=(x_box, y - box_size + 4),
        size=(box_size, box_size),
        fill=color
    ))
    dwg.add(dwg.text(group, insert=(x_text, y), font_size="16px"))
    y += line_spacing

# spacing
y += 40

# ---- SECTION 2: Putative 2ODD Functions (lines) ----
dwg.add(dwg.text("Putative 2ODD function", insert=(30, y), font_size="20px", font_weight="bold"))
y += 40

for func in ordered_functions:
    color = COLORS_2ODD_FUNCTION.get(func, "#cccccc")
    # draw a horizontal line
    dwg.add(dwg.line(
        start=(x_box, y - 6),
        end=(x_box + branch_length, y - 6),
        stroke=color,
        stroke_width=6,  # thick line like tree branches
        stroke_linecap="round"
    ))
    dwg.add(dwg.text(func, insert=(x_text, y), font_size="16px"))
    y += line_spacing

dwg.save()
print("Legend with branch lines created: legend_lines.svg")
# %%
