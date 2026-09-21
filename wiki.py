#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "wikipedia",
#    "truststore"
# ]
# ///

import wikipedia
import truststore
import sys
import requests.exceptions

def lookup_article(query, lang_code):
    """Looks up a Wikipedia article for the given query in the specified language."""
    try:
        print(f"Searching Wikipedia for: {query}...")
        # Set language
        wikipedia.set_lang(lang_code)
        
        # Fetch the page
        page = wikipedia.page(query, auto_suggest=True)
        
        print("\n" + "="*50)
        print(f"Article Title: {page.title}")
        print("="*50)
        print(page.summary)
        print("\n" + "="*50)
        
    except wikipedia.exceptions.PageError:
        print(f"\nError: Wikipedia page for '{query}' not found in language '{lang_code}'.", file=sys.stderr)
    except wikipedia.exceptions.DisambiguationError as e:
        print(f"\nError: Multiple possible matches for '{query}'. Please specify a more precise query.", file=sys.stderr)
        print("Possible options:", e.options[:5], "...", file=sys.stderr)
    except requests.exceptions.SSLError as e:
        print("\n--- SSL/Certificate Error ---", file=sys.stderr)
        print("Could not connect to Wikipedia due to an SSL certificate verification failure.", file=sys.stderr)
        print("This often happens in corporate networks or environments with self-signed certificates.", file=sys.stderr)
        print("Please ensure your system trusts the necessary certificates or configure your environment.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
    except requests.exceptions.RequestException as e:
        print("\n--- Network/Connection Error ---", file=sys.stderr)
        print("Failed to connect to Wikipedia. Check your internet connection, firewall, or proxy settings.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    truststore.inject_into_ssl()
    lang_code = "en" # Default language
    search_args = sys.argv[1:]
    
    # Check for language flag
    if len(search_args) > 0 and search_args[0] in ["-l", "--lang"]:
        lang_code = search_args[1]
        search_args = search_args[2:]
    
    if not search_args:
        print("Usage: python wiki.py [-l <language_code>] <search_term>")
        sys.exit(1)
    
    search_term = " ".join(search_args)
    lookup_article(search_term, lang_code)