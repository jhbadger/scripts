#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "wikipedia-api",
#   "truststore"
# ]
# ///

import sys
import truststore
truststore.inject_into_ssl()

import wikipediaapi

USER_AGENT = "NIHBioinformaticsTool/1.0 (contact: user@nih.gov)"

def lookup_article(query: str, lang_code: str = "en"):
    """Fetches and displays a Wikipedia summary for a given query and language."""
    try:
        print(f"Searching Wikipedia ({lang_code}) for: {query}...")
        
        wiki = wikipediaapi.Wikipedia(
            user_agent=USER_AGENT,
            language=lang_code,
            extract_format=wikipediaapi.ExtractFormat.WIKI
        )
        
        page = wiki.page(query)
        
        if not page.exists():
            print(f"\nError: Article '{query}' does not exist in Wikipedia ({lang_code}).", file=sys.stderr)
            return

        print("\n" + "=" * 50)
        print(f"Article Title: {page.title}")
        print("=" * 50)
        print(page.summary)
        print("\n" + "=" * 50)

    except Exception as e:
        print(f"\nAn error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    lang_code = "en"
    search_args = sys.argv[1:]

    # Correctly grab the language parameter string at index 1
    if len(search_args) >= 2 and search_args[0] in ["-l", "--lang"]:
        lang_code = search_args[1]
        search_args = search_args[2:]

    if not search_args:
        print("Usage: ./wiki.py [-l <language_code>] <search_term>", file=sys.stderr)
        sys.exit(1)

    search_term = " ".join(search_args)
    lookup_article(search_term, lang_code)
