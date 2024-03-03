import matplotlib.pyplot as plt
import matplotlib as mpl


if __name__ == '__main__':
    WGS = True
    mpl.rc('font', family='Times New Roman')

    # Sample data for the pie chart
    labels = ['DeSeq2','ANCOM','ANCOM-BC2', 'LINDA',"LEfSe",'ALDEx']
    sizes = [46.62,13.6,4.56,0.55,32.94,1.73]

    # Create the pie chart
    plt.figure(figsize=(4.5, 4.5))  # Set the figure size if needed
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140,colors=["salmon","mediumaquamarine","green","lightskyblue","gold","grey"])

    # Display the pie chart
    plt.tight_layout()
    plt.show()