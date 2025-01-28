import pandas as pd
import numpy as np
import math

def similar(a, b):
    ### import require packages
    from difflib import SequenceMatcher
    return SequenceMatcher(None, a, b).ratio()

def check_coverage(query, results, similarity=0.75, show_matches=False, show_duplicates=False, show_progress=False):
    '''
    This funcion checks for matches based on title similarity between a reference list of
    article titles and a list of article titles from one or more database searches.

    Parameters:
        - query (list): list of target paper titles to check for.
        - results (list): list of paper titles from database search result.
        - similarity (float): Between 0 and 1, denotes threshold for similarity between titles.
        - show_matches: if True, prints each match along with similarity score.
        - show_duplicates: if True, prints the match every time the target paper is detected
                            more than once in the search (if duplicates remain).
        - show_progress: if True, uses tqdm to display a progress bar.

    Returns:
        List of matches between query list and results list.

    Example:

        matches = check_coverage(query=papers["title"], results=results["title"])

        out:

        There were 52 / 127 unique matches, and 0 duplicates remaining in the search.
        
        
    '''

    ### check threshold is valid
    if not 0 <= similarity <= 1:
                raise ValueError(f"Invalid value for 'similarity'. It should be a float between 0 and 1.")

    print("Checking search results of size", len(results), "for presence of", len(query), "papers...")
    
    n_matches = int(0) # number of unique matches detected
    dupes = int(0) # number of duplicates detected
    matches = [] # empty list for matches

    if show_progress==True:
        from tqdm.autonotebook import tqdm
        for i in tqdm(query): # tqdm displays progress bar
            dup = False
            for j in results:
                if similar(i,j) > similarity:
                    if dup == True:
                        dupes += 1
                        if show_duplicates==True:
                            print("Duplicate detected for:")
                            print(i)
                            print()
        
                    else:
                        dup = True
                        n_matches += 1
                        matches.append(j)
                        if show_matches==True:
                            print("Query:", i)
                            print("Match:", j)
                            print("Similarity:", similar(i,j))
                            print()

    else:
        for i in query:
            dup = False
            for j in results:
                if similar(i,j) > similarity:
                    if dup == True:
                        dupes += 1
                        if show_duplicates==True:
                            print("Duplicate detected for:")
                            print(i)
                            print()
        
                    else:
                        dup = True
                        n_matches += 1
                        matches.append(j)
                        if show_matches==True:
                            print("Query:", i)
                            print("Match:", j)
                            print("Similarity:", similar(i,j))
                            print()
                            
    print("There were", n_matches, "/", len(query), "unique matches, and", dupes, "duplicates remaining in the search.")
    return matches
    