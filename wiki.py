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

def lookup_article(query: str, lang_code: str = "en", show_full: bool = False):
    """Fetches and displays a Wikipedia summary or full article for a given query and language."""
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
        
        if show_full:
            print(page.text)
        else:
            print(page.summary)
            
        print("\n" + "=" * 50)

    except Exception as e:
        print(f"\nAn error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    lang_code = "en"
    show_full = False
    search_term_parts = []
    
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        
        if arg in ["-l", "--lang"]:
            if i + 1 < len(sys.argv):
                lang_code = sys.argv[i+1]
                i += 2
                continue
            else:
                print("Error: -l/--lang requires a language code.", file=sys.stderr)
                sys.exit(1)
        elif arg in ["-f", "--full"]:
            show_full = True
            i += 1
            continue
        else:
            # Assume it's part of the search term
            search_term_parts.append(arg)
            i += 1

    search_term = " ".join(search_term_parts)
    
    if not search_term:
        print("Usage: ./wiki.py [-l <language_code>] [-f|--full] <search_term>", file=sys.stderr)
        sys.exit(1)

    lookup_article(search_term, lang_code, show_full)
