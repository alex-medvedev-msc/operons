def main():
    positive_examples_path = "org_id.csv"

    positive_examples = set()
    with open(positive_examples_path, 'r') as file:
        positive_examples = {line.strip() for line in file.readlines()}

    lacA_examples_path = "lacA.csv"
    lacY_examples_path = "lacY.csv"
    lacZ_examples_path = "lacZ.csv"

    lacA_examples = set()
    lacY_examples = set()
    lacZ_examples = set()

    with open(lacA_examples_path, 'r') as file:
        lacA_examples = {line.strip() for line in file.readlines()}

    with open(lacY_examples_path, 'r') as file:
        lacY_examples = {line.strip() for line in file.readlines()}

    with open(lacZ_examples_path, 'r') as file:
        lacZ_examples = {line.strip() for line in file.readlines()}


    dirty_negative = set.union(lacA_examples, lacY_examples, lacZ_examples)
    clean_negative = dirty_negative.difference(positive_examples)

    with open("negative.csv", 'w') as file:
        file.writelines(str(s)+'\n' for s in clean_negative)


if __name__ == '__main__':
    main()