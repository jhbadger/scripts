# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "fastmcp",
#     "pillow",
# ]
# ///

import json
import os
import base64
from pathlib import Path
from typing import Dict, Any

from fastmcp import FastMCP
from PIL import Image

# Make absolutely sure the root path is fully expanded and resolved immediately
raw_root = os.environ.get("ST_CARDS_DIR", Path.home() / "lib/src/chatbot/cards")
CARDS_ROOT = Path(raw_root).expanduser().resolve()

def get_default_card_path() -> Path:
    """Find a default card if no specific path is given."""
    if CARDS_ROOT.is_dir():
        png_files = list(CARDS_ROOT.glob("*.png"))
        if png_files:
            return png_files[0]
    raise FileNotFoundError(f"No PNG cards found in directory: {CARDS_ROOT}")

def validate_path(path: str) -> Path:
    """Validate and resolve a local card path or partial card name."""
    if not path or not str(path).strip():
        card_path = get_default_card_path()
        return card_path

    given = Path(path).expanduser()

    # If the LLM passes a direct, valid absolute path to a file, use it
    if given.is_absolute() and given.is_file():
        card_path = given.resolve()
    elif CARDS_ROOT.is_dir():
        query = str(path).strip()
        # Search CARDS_ROOT for a partial match
        matches = list(CARDS_ROOT.glob(f"*{query}*.png"))
        if not matches:
            matches = [
                p
                for p in CARDS_ROOT.glob("*.png")
                if query.lower() in p.name.lower()
            ]

        if matches:
            card_path = matches[0].resolve()
        else:
            raise FileNotFoundError(
                f"Could not find any card matching '{query}' in {CARDS_ROOT}"
            )
    else:
        card_path = given.resolve()

    if card_path.suffix.lower() != ".png":
        raise ValueError("The file must be a PNG")

    if not card_path.is_file():
        raise FileNotFoundError(f"Card not found: {card_path}")

    # Ensure it resides in the cards directory
    root_dir = CARDS_ROOT.parent if CARDS_ROOT.is_file() else CARDS_ROOT
    try:
        card_path.relative_to(root_dir)
    except ValueError as exc:
        raise PermissionError(
            f"Card must be inside the configured cards directory: {root_dir}"
        ) from exc

    return card_path

def extract_chara_data(file_path: Path) -> Dict[str, Any]:
    """Extract SillyTavern character data from the PNG tEXt chunk."""
    try:
        with Image.open(file_path) as img:
            img.load()
            if 'chara' in img.info:
                encoded_data = img.info['chara']
                decoded_data = base64.b64decode(encoded_data).decode('utf-8')
                return json.loads(decoded_data)
            else:
                return {"error": "No 'chara' data found in this PNG."}
    except Exception as e:
        return {"error": f"Failed to extract character data: {str(e)}"}

mcp = FastMCP("SillyTavernCards")

@mcp.tool()
def get_card_info(path: str) -> str:
    """Extract and return the text data from a SillyTavern PNG character card. Pass an empty string ("") to use the default card."""
    try:
        card_path = validate_path(path)
        data = extract_chara_data(card_path)
        return json.dumps(data, indent=2)
    except Exception as e:
        return f"Error reading card: {str(e)}"

@mcp.tool()
def get_card_image(path: str) -> str:
    """Gets the character card image. Pass an empty string ("") for the default card."""
    card_path = validate_path(path)
    
    # Read the raw binary data of the image and convert it to a Base64 string
    with open(card_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    
    # Instruct the LLM to output the raw Base64 data as a Markdown image
    return (
        f"Image successfully loaded. To show the image to the user, "
        f"you MUST output exactly this markdown string in your response "
        f"(do not truncate it!):\n\n"
        f"![{card_path.stem}](data:image/png;base64,{encoded_string})"
    )

@mcp.tool()
def list_cards() -> str:
    """List all available PNG character cards in the default directory."""
    if not CARDS_ROOT.is_dir():
        return f"Cards directory not found at {CARDS_ROOT}"
    
    cards = [p.name for p in CARDS_ROOT.glob("*.png")]
    if not cards:
        return f"No PNG cards found in {CARDS_ROOT}"
    
    return f"Found {len(cards)} cards in {CARDS_ROOT}:\n" + "\n".join(cards)

if __name__ == "__main__":
    mcp.run()