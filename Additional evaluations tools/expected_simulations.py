import pandas as pd
import random
def perform_experiment(size, probability):
    # Generate a list of 1s and 0s based on the given probability
    experiment_result = random.choices([1, 0], weights=[probability, 1 - probability], k=size)
    return experiment_result

def load_simulations_params(model,data,size):
    all_info = pd.read_csv(f"rev_fig/consistency_of_model_between_datasets/IBD/{model}.csv",index_col=0)
    if model == "linda" and data == "D1":#or
        num_significant = 0.0
    elif model == "linda-c" and data == "D1":
        num_significant = 0.0
    else:
        num_significant = (all_info[data]).sum()
    probability = num_significant/size
    return probability


if __name__ == '__main__':

    ibd_datasets_names = ["D1", "D2", "D3", "D4", "D5"]
    LIST_MODELS = ["mimic_relative","mimic","linda","linda-c","deseq","deseq-c","ancom","ancom-c","aldex_we","aldex_we-c","aldex_wi","aldex_wi-c","ancom-bc2","ancom-bc2-c","lefse"]
    ITERATIONS =1000
    # LOAD DATA,
    df1 = pd.read_csv("ancom_bc2/data/IBD/tax7_log_subpca.csv", index_col=0)
    df2 = pd.read_csv("ancom_bc2/data/ERP021216/tax7_log_subpca.csv", index_col=0)
    df3 = pd.read_csv("ancom_bc2/data/ok_94/tax7_log_subpca.csv", index_col=0)
    df4 = pd.read_csv("ancom_bc2/data/PRJNA353587/tax7_log_subpca.csv", index_col=0)
    df5 = pd.read_csv("ancom_bc2/data/PRJNA419097/tax7_log_subpca.csv", index_col=0)

    # INTERSECTION
    inter = df1.columns.intersection(df2.columns).intersection(df3.columns).intersection(df4.columns).intersection(
        df5.columns)
    SIZE = len(inter)

    # BUILD EXPECTED PROBABILITY
    for MODEL in LIST_MODELS:
        print(MODEL)
        list_dfs_iters = []
        for iteration in range(ITERATIONS):
            results_df = pd.DataFrame(index=inter,columns=ibd_datasets_names)
            for DATA in ibd_datasets_names:
                probability= load_simulations_params(MODEL,DATA,SIZE)
                results_df[DATA] = perform_experiment(SIZE, probability)
            results_df["Total"] = results_df.sum(axis=1)
            results_df = results_df[results_df["Total"] != 0.0]
            num_s = len(results_df)
            counts = results_df['Total'].value_counts()
            counts = counts.to_frame('Counts')

            # Create a DataFrame with the desired index values
            index_values = [1, 2, 3, 4, 5]
            counts_df = pd.DataFrame(index=index_values)
            counts_df = pd.concat([counts_df, counts], axis=1)
            counts_df = counts_df.fillna(0.0)
            counts_df = (counts_df/num_s)*100
            list_dfs_iters.append(counts_df)

        sum = pd.DataFrame(data=[0.0, 0.0, 0.0, 0.0, 0.0], index=[1, 2, 3, 4, 5], columns=["Counts"])
        for df in list_dfs_iters:
            sum =+df
        sum = sum
        sum.to_csv(f"rev_fig/consistency_of_model_between_datasets/IBD/expected/{MODEL}.csv")
    C=0