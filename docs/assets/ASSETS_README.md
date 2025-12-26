# Assets for GitHub README

This folder contains visual assets for the GitHub repository README.

## Required Assets

### Banner Image
- **File**: `banner.png`
- **Size**: 800 x 200 pixels (recommended)
- **Content**: Project logo with title "AzRM - Azimuthal Reflectivity Method"

### Demo GIF
- **File**: `demo.gif`
- **Size**: 800 x 450 pixels (recommended)
- **Content**: Screen recording of the GUI in action (10-15 seconds loop)
- **Tool**: Use [LICEcap](https://www.cockos.com/licecap/) or [Gifox](https://gifox.io/) to record

### Feature Icons (60x60 pixels each)
- `icon_speed.png` - Lightning bolt or speedometer
- `icon_aniso.png` - Layered earth or wave pattern
- `icon_gui.png` - Browser window or dashboard
- `icon_export.png` - Download arrow or file icon

### Screenshots (800x600 pixels each)
- `screenshot_gather.png` - Seismic gather display (wiggle + image)
- `screenshot_avaz.png` - AVAZ gather display
- `screenshot_avo.png` - AVO/AVAZ curve analysis
- `screenshot_model.png` - Layer editor interface

## Generating Screenshots

### Automated Generation (Recommended)

Run the automated script to generate all screenshots with fixes:

```bash
cd docs
python generate_screenshots.py
```

This script will generate:
- `screenshot_model.png` - Model editor with third layer as OA media
- `screenshot_avo.png` - AVO/AVAZ curves (AVO shows only |Rpp| blue line)
- `screenshot_gather.png` - Seismic gather with fixed polarity at 0° incidence

**Note**: The script includes improved polarity fix for θ=0° using multiple frequency checks.

### Manual Generation

Alternatively, you can capture screenshots manually from the running app:

```bash
cd webapp
streamlit run app.py &
sleep 5

# Use browser developer tools or screenshot tool to capture:
# 1. Navigate to Results tab
# 2. Click "Run Forward Modeling"
# 3. Capture each section
```

## Placeholder Generation

Run this Python script to generate placeholder images:

```python
import matplotlib.pyplot as plt
import numpy as np

# Banner
fig, ax = plt.subplots(figsize=(8, 2))
ax.text(0.5, 0.5, 'AzRM', fontsize=48, fontweight='bold',
        ha='center', va='center', color='#1E88E5')
ax.text(0.5, 0.15, 'Azimuthal Reflectivity Method', fontsize=16,
        ha='center', va='center', color='#666')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')
fig.patch.set_facecolor('#f0f2f6')
plt.savefig('banner.png', dpi=100, bbox_inches='tight',
            facecolor='#f0f2f6', edgecolor='none')
plt.close()
```

## Icon Resources

Free icon sources:
- [Heroicons](https://heroicons.com/)
- [Feather Icons](https://feathericons.com/)
- [Font Awesome](https://fontawesome.com/)
