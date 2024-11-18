import numpy as np
import scipy as sp # may be useful to compute probabilities
import time # may be useful to check the execution time of some function

RedExonPos = np.array([
    [149249757, 149249868], 
    [149256127, 149256423], 
    [149258412, 149258580], 
    [149260048, 149260213], 
    [149261768, 149262007], 
    [149264290, 149264400]  
    ])

GreenExonPos = np.array([
    [149288166, 149288277], 
    [149293258, 149293554], 
    [149295542, 149295710], 
    [149297178, 149297343], 
    [149298898, 149299137], 
    [149301420, 149301530]  
    ])

def loadLastCol(filename):
    global RankMatrix, firstOccurence
    def initialize_firstOccurrence_and_rankMatrix(lastCol):
       
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        n = len(lastCol)
        rankMatrix = {base: np.zeros(n) for base in counts}

        
        for i, char in enumerate(lastCol):
            for base in counts:
                rankMatrix[base][i] = counts[base] + (1 if char == base else 0)
            if char != '$':
                counts[char] += 1
        
        cumulative = 0
        firstOccurrence = {}
        for base in counts:
            firstOccurrence[base] = cumulative
            cumulative += counts[base]
        return firstOccurrence, rankMatrix

    with open(filename) as f:
        lastCol = f.read().replace('\n', '')

    firstOccurence, RankMatrix = initialize_firstOccurrence_and_rankMatrix(lastCol)
    print("First Occurrence and Rank Matrix Created")

    return lastCol  


def MatchReadToLoc(read):
    
    def possible_mismatch_symbols(top, bottom):
       
        return [symbol for symbol in 'ACGT' if RankMatrix[symbol][top - 1] != RankMatrix[symbol][bottom]]

    def calculate_new_bounds(symbol, top, bottom):
       
        new_top = firstOccurence[symbol] + int(RankMatrix[symbol][top]) - (1 if RankMatrix[symbol][top] else 0)
        new_bottom = firstOccurence[symbol] + int(RankMatrix[symbol][bottom]) - (1 if RankMatrix[symbol][bottom] else 0)
        return new_top, new_bottom

    def search_with_mismatch(read, top, bottom, mismatches):
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
        comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join(comp.get(char, 'T') for char in reversed(read))

    positions = []
    positions.extend(search_with_mismatch(read, 0, len(LastCol) - 1, 0))
    positions.extend(search_with_mismatch(generate_reverse_complement(read), 0, len(LastCol) - 1, 0))

    return positions 

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
        if red_counts[i] > 0 and green_counts[i] > 0:  
            overlap_adjustment = (min(red_counts[i], green_counts[i]) / 2)
            adjusted_red[i] -= overlap_adjustment
            adjusted_green[i] -= overlap_adjustment
    return adjusted_red + adjusted_green

def WhichExon(positions):
    
    red_counts = update_exon_counts(positions, RedExonPos)
    green_counts = update_exon_counts(positions, GreenExonPos)
    combined_counts = adjust_counts_for_overlap(red_counts, green_counts)
    
    return np.array(combined_counts)

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
def euclidean_distance(a, b):  
    return np.sqrt(sum((a[i] - b[i]) ** 2 for i in range(len(a))))
def BestMatch(probabilities):
    
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
    LastCol = loadLastCol("../data/chrX_last_col.txt") # loads the last column
    # print("done")
    # print(LastCol.shape())
    RefSeq=open("../data/chrX.fa").read().replace('\n', '').replace(">chrX", '')
    # print("done")
    # print(RefSeq.shape())
    Reads=open("../data/reads").read().strip().split('\n')
    # print("done")
    # print(Reads.shape())
    Map=np.array(open("../data/chrX_map.txt").read().strip().split("\n"), dtype=int)
    print("########### Reads Scanning Started #########")
    # run the functions
    ExonMatchCounts = np.zeros(12) # initialize the counts for exons

    for read in Reads: # update the counts for exons
        positions = MatchReadToLoc(read) # get the list of potential match locations
        # print("aa")
        ExonMatchCounts += WhichExon(positions) # update the counts of exons, if applicable
    # print("done")
    ListProb = ComputeProb(ExonMatchCounts) # compute probabilities of each of the four configurations
    MostLikely = BestMatch(ListProb) # find the most likely configuration
    print("Configuration %d is the best match"%MostLikely)
