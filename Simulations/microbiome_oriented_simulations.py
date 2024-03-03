import pandas as pd
import numpy as np
import random


def increase_in_p(df, bact_name, p):
    list_indexes = [i for i in df.columns]
    numerical_columns = [word for word in list_indexes if word != "taxonomy"]
    df[numerical_columns] = df[numerical_columns].astype(float)
    print(df[numerical_columns].sum().sum())
    for t, df_t in df.groupby("taxonomy"):
        if t == bact_name:
            # Update numerical columns with increased values
            df.loc[df_t.index, numerical_columns] *= (100 + p) / 100
    print(df[numerical_columns].sum().sum())
    return df


if __name__ == '__main__':
    COMMON = True
    LIST_NAMES = ["IBD"]  # ["ERP020401","IBD","PRJNA353587"]
    for NAME in LIST_NAMES:
        df = pd.read_csv(f"Data_for_git_process/{NAME}/for_preprocess.csv", index_col=0)
        processed = pd.read_csv(f"corr_results/deseq/data_for_deseq_{NAME}.csv", index_col=0)
        n = int(len(processed.index) * 0.2)
        to_choose = max(20, n)

        list_indexes = [i for i in df.index]
        list_indexes_without_taxonomy = [word for word in list_indexes if word != "taxonomy"]

        # DEFINE A RADOM TAG
        num_samples = len(df) - 1
        half_samples = num_samples // 2
        half_samples_ = num_samples - half_samples

        # Generate a list of zeros and ones
        tags = np.concatenate((np.zeros(half_samples), np.ones(half_samples_)))

        # Shuffle the tags
        np.random.shuffle(tags)
        tag = pd.DataFrame(index=list_indexes_without_taxonomy, columns=["Tag"], data=tags)

        # CHOOSE 10 BACTS FOR NEGATIVE
        if COMMON:
            num_non_zeros = (processed != 0.0).sum()
            num_non_zeros = num_non_zeros.sort_values(ascending=False)
            if NAME == "PRJNA353587":
                del num_non_zeros['Unassigned']
            num_non_zeros = num_non_zeros.head(to_choose)
            my_list = list(num_non_zeros.index)
        else:
            processed = processed.sum()
            processed_non_zero = processed[processed.values != 0.0]
            my_list = processed.index
            my_list = list(set(my_list))
        # Choose 10 entries randomly from the list
        random_entries_1 = random.sample(my_list, k=10)

        # Remove the selected entries from the original list
        remaining_entries = [entry for entry in my_list if entry not in random_entries_1]

        # Choose 10 more entries randomly from the remaining list
        random_entries_2 = random.sample(remaining_entries, k=10)
        print(random_entries_1)
        print(random_entries_2)

        # DIVIDE SAMPLES ACCORDING TO TAG
        df = df.T
        tag0 = tag[tag["Tag"] == 0]
        tag1 = tag[tag["Tag"] == 1]
        df0 = df[list(tag0.index) + ["taxonomy"]]
        df1 = df[list(tag1.index) + ["taxonomy"]]

        # INCREASE AMOUNT
        for b in random_entries_1:
            df0 = increase_in_p(df0, b, 20)
        for b in random_entries_2:
            df1 = increase_in_p(df1, b, 20)

    to_save = pd.concat([df0.T, df1.T])
    to_save.to_csv(f"corr_results/sim_{NAME}/common/for_preprocess.csv")
    tag.to_csv(f"corr_results/sim_{NAME}/common/tag.csv")
