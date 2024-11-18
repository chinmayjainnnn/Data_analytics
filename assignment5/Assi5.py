import numpy as np
import scipy as sp # may be useful to compute probabilities
import time # may be useful to check the execution time of some function





RedExonPos = np.array([
    [149249757, 149249868], # R1
    [149256127, 149256423], # R2
    [149258412, 149258580], # R3
    [149260048, 149260213], # R4
    [149261768, 149262007], # R5
    [149264290, 149264400]  # R6
    ])

GreenExonPos = np.array([
    [149288166, 149288277], # G1
    [149293258, 149293554], # G2
    [149295542, 149295710], # G3
    [149297178, 149297343], # G4
    [149298898, 149299137], # G5
    [149301420, 149301530]  # G6
    ])

def loadLastCol(filename):
    """
    Input: Path of the file corresponding to the last column (BWT).
    Output: The last column (BWT) in string format.
    """
    global RankMatrix, firstOccurence

    def initialize_firstOccurrence_and_rankMatrix(lastCol):
        # Initialize counters and rank matrix
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        n = len(lastCol)
        rankMatrix = {base: np.zeros(n) for base in counts}

        # Populate Rank Matrix based on character frequencies
        for i, char in enumerate(lastCol):
            for base in counts:
                rankMatrix[base][i] = counts[base] + (1 if char == base else 0)
            if char != '$':
                counts[char] += 1

        # Calculate first occurrence for each character
        cumulative = 0
        firstOccurrence = {}
        for base in counts:
            firstOccurrence[base] = cumulative
            cumulative += counts[base]

        return firstOccurrence, rankMatrix

    # Load the last column from the file and initialize structures
    with open(filename) as f:
        lastCol = f.read().replace('\n', '')

    firstOccurence, RankMatrix = initialize_firstOccurrence_and_rankMatrix(lastCol)
    print("First Occurrence and Rank Matrix Created")

    return lastCol  # Return last column as a string





# def MatchReadToLoc(read):
#     """
#     Input: A read (string).
#     Output: List of potential locations at which the read may match the reference sequence.
#     """
#     positions = []

#     def possible_mismatch_symbols(top, bottom):
#         return [symbol for symbol in 'ACGT' if RankMatrix[symbol][top - 1] != RankMatrix[symbol][bottom]]

#     def search_with_mismatch(read, top, bottom, mismatches):
#         if top > bottom:
#             return []

#         while read:
#             symbol = read[-1] if read[-1] != 'N' else 'A'
#             read = read[:-1]

#             if RankMatrix[symbol][top - 1] == RankMatrix[symbol][bottom] and top > 0:
#                 if mismatches < 2:
#                     mismatches += 1
#                     match_positions = []
#                     for alt_symbol in possible_mismatch_symbols(top, bottom):
#                         new_top = firstOccurence[alt_symbol] + int(RankMatrix[alt_symbol][top]) - (1 if RankMatrix[alt_symbol][top] else 0)
#                         new_bottom = firstOccurence[alt_symbol] + int(RankMatrix[alt_symbol][bottom]) - (1 if RankMatrix[alt_symbol][bottom] else 0)
#                         match_positions.extend(search_with_mismatch(read, new_top, new_bottom, mismatches))
#                     return match_positions
#                 else:
#                     return []
#             else:
#                 top = firstOccurence[symbol] + int(RankMatrix[symbol][top]) - (1 if RankMatrix[symbol][top] else 0)
#                 bottom = firstOccurence[symbol] + int(RankMatrix[symbol][bottom]) - (1 if RankMatrix[symbol][bottom] else 0)
#         return [int(Map[i]) for i in range(top + 1, bottom + 1)]

#     # Generate reverse complement of the read
#     comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
#     rev_comp_read = ''.join(comp.get(char, 'T') for char in reversed(read))

#     # Perform search for both the read and its reverse complement
#     positions.extend(search_with_mismatch(read, 0, len(LastCol) - 1, 0))
#     positions.extend(search_with_mismatch(rev_comp_read, 0, len(LastCol) - 1, 0))
#     # print("ii")
#     return positions  # List of potential match locations

def MatchReadToLoc(read):
    """
    Input: A read (string).
    Output: List of potential locations at which the read may match the reference sequence.
    """
    def possible_mismatch_symbols(top, bottom):
        """Returns possible mismatch symbols between top and bottom positions in RankMatrix."""
        return [symbol for symbol in 'ACGT' if RankMatrix[symbol][top - 1] != RankMatrix[symbol][bottom]]

    def calculate_new_bounds(symbol, top, bottom):
        """Calculates new top and bottom bounds for the given symbol."""
        new_top = firstOccurence[symbol] + int(RankMatrix[symbol][top]) - (1 if RankMatrix[symbol][top] else 0)
        new_bottom = firstOccurence[symbol] + int(RankMatrix[symbol][bottom]) - (1 if RankMatrix[symbol][bottom] else 0)
        return new_top, new_bottom

    def search_with_mismatch(read, top, bottom, mismatches):
        """Performs search allowing up to 1 mismatch in the read."""
        if top > bottom:
            return []

        positions = []
        while read:
            symbol = read[-1] if read[-1] != 'N' else 'A'
            read = read[:-1]

            if RankMatrix[symbol][top - 1] == RankMatrix[symbol][bottom] and top > 0:
                if mismatches < 2:
                    mismatches += 1
                    for alt_symbol in possible_mismatch_symbols(top, bottom):
                        new_top, new_bottom = calculate_new_bounds(alt_symbol, top, bottom)
                        positions.extend(search_with_mismatch(read, new_top, new_bottom, mismatches))
                    return positions
                else:
                    return []
            else:
                top, bottom = calculate_new_bounds(symbol, top, bottom)
        
        positions.extend([int(Map[i]) for i in range(top + 1, bottom + 1)])
        return positions

    def generate_reverse_complement(read):
        """Generates the reverse complement of a DNA sequence."""
        comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join(comp.get(char, 'T') for char in reversed(read))

    # Perform search for both the read and its reverse complement
    positions = []
    positions.extend(search_with_mismatch(read, 0, len(LastCol) - 1, 0))
    positions.extend(search_with_mismatch(generate_reverse_complement(read), 0, len(LastCol) - 1, 0))

    return positions  # List of potential match locations

# def WhichExon(positions):
#     """
#     Input: list of potential locations at which the read may match the reference sequence.
#     Output: Update(increment) to the counts of the 12 exons
#     """
    
#     resR = [0]*6
#     resG = [0]*6
#     flag= [0]*6
    
#     for pos in positions:             
#         for i in range(6):
#             if pos >= RedExonPos[i][0] and pos <= RedExonPos[i][1]:   #update if any Red Exon location matches
#                 resR[i]+=1
#                 flag[i]+=1
#             if pos >= GreenExonPos[i][0] and pos <= GreenExonPos[i][1]: #update if any Red Exon location matche
#                 resG[i]+=1
#                 flag[i]+=1
    
#     for i in range(6):                 #used to update values that are updated for both green and red for a particular exon 
#         if flag[i]>1 :                 # for both red and green correspondingly
#             resR[i] -= ((flag[i]//2)*1/2)
#             resG[i] -= ((flag[i]//2)*1/2)

#     resR = resR + resG
    
#     return np.array(resR)
def update_exon_counts(positions, exon_positions):

    counts = [0] * len(exon_positions)
    for pos in positions:
        for i, (start, end) in enumerate(exon_positions):
            if start <= pos <= end:
                counts[i] += 1
    return counts

def adjust_counts_for_overlap(red_counts, green_counts):
  
    adjusted_red = red_counts[:]
    adjusted_green = green_counts[:]
    for i in range(len(red_counts)):
        if red_counts[i] > 0 and green_counts[i] > 0:  # if both red and green have counts
            overlap_adjustment = (min(red_counts[i], green_counts[i]) / 2)
            adjusted_red[i] -= overlap_adjustment
            adjusted_green[i] -= overlap_adjustment
    return adjusted_red + adjusted_green

def WhichExon(positions):
    
    red_counts = update_exon_counts(positions, RedExonPos)
    green_counts = update_exon_counts(positions, GreenExonPos)
    combined_counts = adjust_counts_for_overlap(red_counts, green_counts)
    
    return np.array(combined_counts)



# def ComputeProb(ExonMatchCounts):
#     """
#     Input: The counts for each exon
#     Output: Probabilities of each of the four configurations (a list of four real numbers)
#     """
#     print("Exons Counts: ",ExonMatchCounts)
    
    
#     # function body - Begins
#     P0 = ExonMatchCounts[1]/ExonMatchCounts[7]
#     P1 = ExonMatchCounts[2]/ExonMatchCounts[8]
#     P2 = ExonMatchCounts[3]/ExonMatchCounts[9]
#     P3 = ExonMatchCounts[4]/ExonMatchCounts[10]
#     # function body - ends
#     print("Probablity of Given Exons is as follows: ",[P0, P1, P2, P3])
#     return [P0, P1, P2, P3]


def calculate_probability(count, total):
    return count / total if total != 0 else 0

def extract_counts(exon_match_counts):
    return [
        (exon_match_counts[1], exon_match_counts[7]),
        (exon_match_counts[2], exon_match_counts[8]),
        (exon_match_counts[3], exon_match_counts[9]),
        (exon_match_counts[4], exon_match_counts[10])
    ]

def ComputeProb(exon_match_counts):

    counts = extract_counts(exon_match_counts)
    probabilities = [calculate_probability(count, total) for count, total in counts]

    print("Exon Counts:", exon_match_counts)
    print("Probabilities of Given Exons:", probabilities)
    
    return probabilities

def BestMatch(probabilities):
    
    def euclidean_distance(a, b):
       
        return np.sqrt(sum((a[i] - b[i]) ** 2 for i in range(len(a))))

    def find_best_configuration(probabilities, configurations):
        
        min_distance = float('inf')
        best_config_index = -1
        for i, config in enumerate(configurations):
            distance = euclidean_distance(probabilities, config)
            if distance < min_distance:
                min_distance = distance
                best_config_index = i
        return best_config_index

    
    configurations = [
        [0.5, 0.5, 0.5, 0.5],
        [1, 1, 0, 0],
        [0.33, 0.33, 1, 1],
        [0.33, 0.33, 0.33, 1]
    ]

    return find_best_configuration(probabilities, configurations)


if __name__ == "__main__":
    # load all the data files
    LastCol = loadLastCol("../data/chrX_last_col.txt") # loads the last column
    print("done")
    RefSeq=open("../data/chrX.fa").read().replace('\n', '').replace(">chrX", '')
    print("done")
    Reads=open("../data/reads").read().strip().split('\n')
    print("done")
    # Map = loadMapToRefSeq("../data/chrX_map.txt") # loads the mapping to the reference sequence
    Map=np.array(open("../data/chrX_map.txt").read().strip().split("\n"), dtype=int)
    print("########### Map Loaded #########")
    print("########### Reads Scanning Started #########")
    # run the functions
    ExonMatchCounts = np.zeros(12) # initialize the counts for exons

    for read in Reads: # update the counts for exons
        positions = MatchReadToLoc(read) # get the list of potential match locations
        print("aa")
        ExonMatchCounts += WhichExon(positions) # update the counts of exons, if applicable
    print("done")
    ListProb = ComputeProb(ExonMatchCounts) # compute probabilities of each of the four configurations
    MostLikely = BestMatch(ListProb) # find the most likely configuration
    print("Configuration %d is the best match"%MostLikely)
