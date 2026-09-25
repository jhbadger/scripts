#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "fastmcp>=2.0",
#   "pillow>=10.0",
# ]
# ///

from __future__ import annotations

import base64
import json
from pathlib import Path
from typing import Any

from fastmcp import FastMCP
from PIL import Image

mcp = FastMCP("SillyTavern Character Cards")

# Optional safety restriction.
# Set this to the directory containing your cards.
# Leave as None to allow any local PNG path.
CARDS_ROOT = Path.home() / "lib/src/chatbot/cards"


def validate_path(path: str) -> Path:
    """Validate and resolve a local card path."""
    card_path = Path(path).expanduser().resolve()

    if card_path.suffix.lower() != ".png":
        raise ValueError("The file must be a PNG")

    if not card_path.is_file():
        raise FileNotFoundError(f"Card not found: {card_path}")

    if CARDS_ROOT is not None:
        root = CARDS_ROOT.expanduser().resolve()

        try:
            card_path.relative_to(root)
        except ValueError as exc:
            raise PermissionError(
                f"Card must be inside the configured cards directory: {root}"
            ) from exc

    return card_path


def read_card(path: str) -> dict[str, Any]:
    """Extract the embedded SillyTavern character JSON."""
    card_path = validate_path(path)

    with Image.open(card_path) as image:
        encoded = image.info.get("chara")

    if not encoded:
        raise ValueError(
            "This PNG does not contain SillyTavern 'chara' metadata"
        )

    try:
        if isinstance(encoded, bytes):
            encoded = encoded.decode("ascii")

        decoded = base64.b64decode(encoded)
        card = json.loads(decoded.decode("utf-8"))
    except Exception as exc:
        raise ValueError(f"Could not decode character card: {exc}") from exc

    if not isinstance(card, dict):
        raise ValueError("Character-card data is not a JSON object")

    return card


def card_data(card: dict[str, Any]) -> dict[str, Any]:
    """Support both V2 cards and older flat card formats."""
    data = card.get("data")

    if isinstance(data, dict):
        return data

    return card


@mcp.tool
def inspect_card(path: str) -> dict[str, Any]:
    """
    Inspect a SillyTavern PNG character card.

    Returns the character information without the image itself.
    """
    card = read_card(path)
    data = card_data(card)

    return {
        "spec": card.get("spec"),
        "spec_version": card.get("spec_version"),
        "name": data.get("name"),
        "description": data.get("description"),
        "personality": data.get("personality"),
        "scenario": data.get("scenario"),
        "first_mes": data.get("first_mes"),
        "alternate_greetings": data.get("alternate_greetings", []),
        "mes_example": data.get("mes_example"),
        "system_prompt": data.get("system_prompt"),
        "post_history_instructions": data.get(
            "post_history_instructions"
        ),
        "tags": data.get("tags", []),
        "creator": data.get("creator"),
        "creator_notes": data.get("creator_notes"),
        "character_version": data.get("character_version"),
        "has_character_book": bool(data.get("character_book")),
        "extensions": data.get("extensions", {}),
    }


@mcp.tool
def load_card(path: str) -> str:
    """
    Convert a SillyTavern card into a prompt suitable for a roleplay session.
    """
    card = read_card(path)
    data = card_data(card)

    sections: list[str] = []

    def add_section(title: str, value: Any) -> None:
        if value:
            sections.append(f"## {title}\n{value}")

    add_section("Character Name", data.get("name"))
    add_section("Description", data.get("description"))
    add_section("Personality", data.get("personality"))
    add_section("Scenario", data.get("scenario"))
    add_section("Example Dialogue", data.get("mes_example"))
    add_section("System Prompt", data.get("system_prompt"))
    add_section(
        "Post-History Instructions",
        data.get("post_history_instructions"),
    )

    character_book = data.get("character_book")

    if isinstance(character_book, dict):
        entries = character_book.get("entries", [])

        if entries:
            lore = []

            for entry in entries:
                if not isinstance(entry, dict):
                    continue

                content = entry.get("content")
                keys = entry.get("keys", [])

                if content:
                    lore.append(
                        f"Keys: {', '.join(map(str, keys))}\n"
                        f"Content: {content}"
                    )

            if lore:
                add_section("Character Lorebook", "\n\n".join(lore))

    return "\n\n".join(sections)


@mcp.tool
def get_greeting(
    path: str,
    alternate_index: int | None = None,
) -> str:
    """
    Return the card's first greeting or an alternate greeting.
    """
    card = read_card(path)
    data = card_data(card)

    if alternate_index is None:
        return str(data.get("first_mes", ""))

    greetings = data.get("alternate_greetings", [])

    if not isinstance(greetings, list):
        raise ValueError("alternate_greetings is not a list")

    if not 0 <= alternate_index < len(greetings):
        raise IndexError(
            f"alternate_index must be between 0 and {len(greetings) - 1}"
        )

    return str(greetings[alternate_index])


@mcp.tool
def list_cards(directory: str) -> list[str]:
    """List PNG files in a local character-card directory."""
    folder = Path(directory).expanduser().resolve()

    if not folder.is_dir():
        raise NotADirectoryError(f"Not a directory: {folder}")

    return sorted(
        str(path)
        for path in folder.glob("*.png")
        if path.is_file()
    )


if __name__ == "__main__":
    mcp.run()
