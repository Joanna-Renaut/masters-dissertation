import codecs
import csv
import time
import xml.etree.ElementTree as etree

# useful information and code
# https://docs.drugbank.com/xml/#introduction
# https://github.com/jeffheaton/article-code/blob/master/python/wikipedia/wiki-basic-stream.py


xml_file = "drugbank_full.xml"
csv_file = "../database-csv/all_drugbank.csv"
encoding = "utf-8"


# # # time string


def time_string(sec_elapsed):
    h = int(sec_elapsed / (60 * 60))
    m = int((sec_elapsed % (60 * 60)) / 60)
    s = sec_elapsed % 60
    return "{}:{:>02}:{:>05.2f}".format(h, m, s)


# # strip the tag name


def strip_tag_name(t):
    t = elem.tag
    idx = k = t.rfind("}")
    if idx != -1:
        t = t[idx + 1 :]
    return t


def is_child_of(tag_name, path):
    if len(path) < 1:  # if path is less than 1, stop
        return False

    parent = path[-1]  # the -1 position in the path is the parent
    # print(parent)
    return tag_name == parent


drugs = etree.iterparse(xml_file, events=("start", "end"))

drug_count = 0
name = None
start_time = time.time()
path = []  # create an empty list for the paths to go to

# # set up for writing to csv

with codecs.open(csv_file, "w", encoding) as drugbank_drugs:
    drugs_writer = csv.writer(drugbank_drugs, quoting=csv.QUOTE_MINIMAL)
    drugs_writer.writerow(
        [
            "name",
            "description",
            "state",
            "synonyms",
            "indications",
            "dosages",
            "drug interactions",
            "pathways",
            "target",
        ]
    )

    for event, elem in drugs:
        tag_name = strip_tag_name(elem.tag)

        # # if there is a start event (e.g. <drug>) append the path list up to the tag name
        # e.g. tag_name = name # path = ['drugbank', 'drug', 'name']
        # e.g. tag_name = description # path = ['drugbank', 'drug', 'description']

        if event == "start":
            path.append(tag_name)
        else:
            path.pop()  # removes the tag_name at the end e.g. <drug> <name> to just <drug>
        # print(path)

        # # check the start event, reset all the values and then check the end event - if wanted, write element
        if event == "start":
            if tag_name == "drug":
                drug_name = ""
                description = ""
                state = ""
                synonyms = []
                indication = ""
                dosages = []
                drug_interactions = []
                pathways = []
                targets = []

        else:  # event == 'end'
            if tag_name == "name":
                if is_child_of(
                    "drug", path
                ):  # if it is a child of the tag name (e.g. name within drug) then add text
                    drug_name = elem.text
                    # print(elem.text)
                elif is_child_of("drug-interaction", path):
                    drug_interaction_name = elem.text
                elif is_child_of("pathway", path):
                    pathway_name = elem.text
            elif tag_name == "description":
                if is_child_of("drug", path):
                    description = elem.text
                    # print(elem.text)
            elif tag_name == "state":
                state = elem.text
            elif tag_name == "synonym":
                if is_child_of("synonyms", path):
                    synonym = elem.text
            elif tag_name == "synonyms":
                synonyms.append(synonym)
                # print(synonyms)
            elif tag_name == "indication":
                indication = elem.text
                # print(elem.text)
            elif tag_name == "form":
                form = elem.text
                # print(elem.text)
            elif tag_name == "route":
                route = elem.text
                # print(elem.text)
            elif tag_name == "strength":
                strength = elem.text
                # print(elem.text)
            elif tag_name == "dosage":
                dosages.append(form)
                dosages.append(route)
                dosages.append(strength)
            elif tag_name == "drug-interactions":
                drug_interactions.append(drug_interaction_name)
                # print(drug_interactions)
            elif tag_name == "pathways":
                pathways.append(pathway_name)
                # print(pathways)
            elif tag_name == "gene-name":
                gene_name = elem.text
                # print(elem.text)
            elif tag_name == "polypeptide":
                targets.append(gene_name)
            elif tag_name == "drug":
                drug_count += 1

                # print(drug_name)
                # print(targets)
                # print(description)

                # # output a csv row for each gene name / target
                # # duplicates drugs, but will be easier to compare SSL csv etc
                for target in targets:
                    drugs_writer.writerow(
                        [
                            drug_name,
                            description,
                            state,
                            synonyms,
                            indication,
                            dosages,
                            drug_interactions,
                            pathways,
                            target,
                        ]
                    )

            elem.clear()  # clear = nice performance and lower memory usage

elapsed_time = time.time() - start_time

print(f"Total drugs = {drug_count}")
print(f"Total time = {elapsed_time}")
