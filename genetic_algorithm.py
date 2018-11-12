
#!/usr/bin/env python
import random
import datetime

geneSet = "AGTC"
target = "CCGCCCGGGGTTGAACGATGCATTAAGCGTATATAGACGGTTGTGCATGTCCACTGGAACAAGGCATCCCCACCTACACGACCCTTGCAAGGAACCTAGA"

#generate a parental string of letters to mutate
def generate_parent(length):
    genes = []
    while len(genes) < length:
        sampleSize = min(length - len(genes), len(geneSet))
        genes.extend(random.sample(geneSet, sampleSize))
    return ''.join(genes)

print("something is running")

# feedback to the engine to guide ot toward a solution for a given sequence called guess
# Total number of letters in the guess sequence that match the letter in the same position of the target

def get_fitness(guess):
    return sum(1 for expected, actual in zip(target, guess)
               if expected == actual)

# Convert the parent string to an array, then replaces one letter with another randomly selected
def mutate(parent):
    index = random.randrange(0, len(parent))
    childGenes = list(parent)
    newGene, alternate = random.sample(geneSet, 2)  # Extract two random letters from geneSet and store in two different
                                                    # variables: newGene and alternate (is a list with two values)

    # split a command that is originally in one line in multiple lines (\)
    childGenes[index] = alternate \
        if newGene == childGenes[index] \
        else newGene
    return ''.join(childGenes)


def display(guess):
    timeDiff = datetime.datetime.now() - startTime
    fitness = get_fitness(guess)
    print("{0}\t{1}\t{2}".format(guess, fitness, str(timeDiff)))

random.seed()
startTime = datetime.datetime.now()
bestParent = generate_parent(len(target))
bestFitness = get_fitness(bestParent)
display(bestParent)

while True:
    child = mutate(bestParent)
    childFitness = get_fitness(child)

    if bestFitness >= childFitness:
        continue
    display(child)
    if childFitness >= len(bestParent):
        break
    bestFitness = childFitness
    bestParent = child
