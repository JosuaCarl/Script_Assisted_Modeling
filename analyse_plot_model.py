import sys
import os
import libsbml
import matplotlib.pyplot as plt
import numpy as np


# Usage: analyse_plot_model.py <path_1> <path_2> ... <label_1> <label_2> ...
# Number of labels must be equal to number of paths !
def main(args):
    # console access
    if len(args) % 2 != 1:
        print(main.__doc__)
        sys.exit(1)

    hp = int((len(args) - 1) / 2) + 1
    model_paths = args[1:hp]
    program_labels = args[hp:]

    reader = libsbml.SBMLReader()

    for mp in model_paths:
        if not os.path.exists(mp):
            print("[Error] %s : No such file." % mp)
            sys.exit(1)

    # import model in sbml/xml format'
    # Convert XML to JSON format (for  Escher)
    # json_path = outpath + program_name + "_" + program_version + ".json"
    # cobra.io.save_json_model(model, filename=json_path)

    # Escher Build for Visualization
    # builder = escher.Builder(model_json=json_path)

    # import docs
    docs = []
    for mp in model_paths:
        docs.append(reader.readSBML(mp))
    print(docs)

    # Extract model info
    models = []
    for doc in docs:
        models.append(doc.getModel())

    # Extract reactions, metabolites & genes
    reactions = []
    metabolites = []
    genes = []
    for m in models:
        reactions.append(m.getNumReactions())
        metabolites.append(m.getNumSpecies())
        genes.append(m.getPlugin('fbc').getNumGeneProducts())

    # Process Data for usage in plot
    models_data = []
    for i in range(0, len(reactions)):
        models_data.append([reactions[i], metabolites[i], genes[i]])
    labels = ["reactions", "metabolites", "genes"]
    x = np.arange(len(labels))  # the label locations
    num_data = len(models_data)
    width = 0.7 / num_data  # the width of the bars
    place_bars = width / num_data

    # Plotting
    rects = []
    fig, ax = plt.subplots()
    for j in range(0, len(models_data)):
        x_position = x - width + place_bars * num_data * (j + 1)
        rects.append(ax.bar(x_position, models_data[j], width, label=program_labels[j]))

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_title('Models in comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    for rect in rects:
        ax.bar_label(rect, padding=3)

    fig.tight_layout()

    plt.show()


if __name__ == '__main__':
    main(sys.argv)
