import warnings

import numpy as np
import pandas as pd
import scipy.io
from sklearn import svm
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.model_selection import train_test_split
from tqdm import tqdm

from graces import GRACES

warnings.filterwarnings("ignore")

# grid search range on key hyperparameters
dropout_prob = [0.1, 0.5, 0.75]
f_correct = [0, 0.1, 0.5, 0.9]

# declare empty dataframe
results = {
    'number_of_features': [], 'iteration': [], 'feature_selection_indices': [],
    'feature_selection_protein': [], 'roc_auc_score': [], 'average_roc_auc_score': [],
    'average_testing_accuracy': [], 'max_testing_accuracy': [], 'max_testing_accuracy_rs': [],
    'most_mismatch_data_sample_test': [], 'average_full_dataset_accuracy': [],
    'max_full_dataset_accuracy': [], 'max_full_dataset_accuracy_rs': [], 'most_mismatch_data_sample_full': []
}


def main(name, n_features, n_iters, n_repeats):
    # np.random.seed(0)  # for reproducibility
    data = scipy.io.loadmat('annotation/' + name)
    x = data['X'].astype(float)
    y = np.int64(data['Y'])
    y = y.reshape(-1)
    p_names = data["all_proteins"]
    b_names = data["bacteria_name"]
    # auc_test = np.zeros(n_iters)
    seeds = np.random.choice(range(100), n_iters, replace=False)  # for reproducibility
    for iter in tqdm(range(n_iters), desc="Iteration Variation"):
        x_train, x_test, y_train, y_test = train_test_split(x, y, train_size=0.7, random_state=seeds[iter], stratify=y)
        x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, train_size=2 / 7, random_state=seeds[iter],
                                                          stratify=y_train)
        auc_grid = np.zeros((len(dropout_prob), len(f_correct)))
        loss_grid = np.zeros((len(dropout_prob), len(f_correct)))
        for i in range(len(dropout_prob)):
            for j in range(len(f_correct)):
                for r in range(n_repeats):
                    slc_g = GRACES(n_features=n_features, dropout_prob=dropout_prob[i], f_correct=f_correct[j])
                    selection_g = slc_g.select(x_train, y_train)
                    x_train_red_g = x_train[:, selection_g]
                    x_val_red = x_val[:, selection_g]
                    clf_g = svm.SVC(probability=True)
                    clf_g.fit(x_train_red_g, y_train)
                    y_val_pred = clf_g.predict_proba(x_val_red)
                    auc_grid[i, j] += roc_auc_score(y_val, y_val_pred[:, 1])
                    loss_grid[i, j] += -np.sum(y_val * np.log(y_val_pred[:, 1]))
        index_i, index_j = np.where(auc_grid == np.max(auc_grid))
        best_index = np.argmin(loss_grid[index_i, index_j])  # break tie based on cross-entropy loss
        best_prob, best_f_correct = dropout_prob[int(index_i[best_index])], f_correct[int(index_j[best_index])]
        auc_score_arr = []
        selection_arr = []
        for r in range(n_repeats):
            slc = GRACES(n_features=n_features, dropout_prob=best_prob, f_correct=best_f_correct)
            selection = slc.select(x_train, y_train)
            selection_arr.append(selection)
            x_train_red = x_train[:, selection]
            x_test_red = x_test[:, selection]
            clf = svm.SVC(probability=True)
            clf.fit(x_train_red, y_train)
            y_test_pred = clf.predict_proba(x_test_red)
            auc_score_arr.append(roc_auc_score(y_test, y_test_pred[:, 1]))
        append_to_results(results, 'number_of_features', n_features)
        append_to_results(results, 'iteration', iter)
        append_to_results(results, 'feature_selection_indices', selection_arr[np.argmax(auc_score_arr)])
        append_to_results(results, 'feature_selection_protein', p_names[selection_arr[np.argmax(auc_score_arr)]])
        append_to_results(results, 'roc_auc_score', np.max(auc_score_arr))
        append_to_results(results, 'average_roc_auc_score', np.average(auc_score_arr))
        test_on_feature_set(selection_arr[np.argmax(auc_score_arr)], x=x, y=y, b_names=b_names)


def test_on_feature_set(feature_selection, name=None, x=None, y=None, b_names=None):
    if name:
        data = scipy.io.loadmat('annotation/' + name)
        x = data['X'].astype(float)
        y = np.int64(data['Y'])
        y = y.reshape(-1)
        b_names = data["bacteria_name"]
    test_accuracies = []
    full_accuracies = []
    mismatch_full_data_samples = []
    mismatch_test_data_samples = []
    for i in range(10000):
        x_selected = x[:, feature_selection]
        x_train, x_test, y_train, y_test = train_test_split(x_selected, y, train_size=0.7, random_state=i,
                                                            stratify=y)
        b_train, b_test, _, _ = train_test_split(b_names, y, train_size=0.7, random_state=i, stratify=y)
        clf = svm.SVC()
        clf.fit(x_train, y_train)
        y_test_pred = clf.predict(x_selected)
        mismatch_full_data_samples.extend(b_names[find_mismatch_indices(y_test_pred, y)])
        accuracy = accuracy_score(y, y_test_pred)
        full_accuracies.append(accuracy)
        y_test_pred = clf.predict(x_test)
        accuracy = accuracy_score(y_test, y_test_pred)
        test_accuracies.append(accuracy)
        mismatch_test_data_samples.extend(b_test[find_mismatch_indices(y_test_pred, y_test)])

    append_to_results(results, 'average_testing_accuracy', np.average(test_accuracies))
    append_to_results(results, 'max_testing_accuracy', np.max(test_accuracies))
    append_to_results(results, 'max_testing_accuracy_rs', np.argmax(test_accuracies))
    append_to_results(results, 'most_mismatch_data_sample_test', most_common_element(mismatch_test_data_samples))
    append_to_results(results, 'average_full_dataset_accuracy', np.average(full_accuracies))
    append_to_results(results, 'max_full_dataset_accuracy', np.max(full_accuracies))
    append_to_results(results, 'max_full_dataset_accuracy_rs', np.argmax(full_accuracies))
    append_to_results(results, 'most_mismatch_data_sample_full', most_common_element(mismatch_full_data_samples))


def most_common_element(lst):
    """
    Finds and returns the most common element in the list using NumPy.

    Parameters:
        lst (list): The input list.

    Returns:
        The most common element in the list.
    """
    if not lst:
        return None

    unique_elements, counts = np.unique(lst, return_counts=True)
    max_count_index = np.argmax(counts)

    return unique_elements[max_count_index]


def find_mismatch_indices(list1, list2):
    if len(list1) != len(list2):
        raise ValueError("Input lists must have the same length")

    mismatch_indices = [index for index, (elem1, elem2) in enumerate(zip(list1, list2)) if elem1 != elem2]
    return mismatch_indices


def append_to_results(dict, column_name, value):
    dict[column_name].append(value)


if __name__ == "__main__":
    name = 'pseudomonas_data'
    max_features = 70
    n_iters = 30
    n_repeats = 3

    for i in tqdm(range(1, max_features + 1), desc="Feature Variation"):
        main(name=name, n_features=i, n_iters=n_iters, n_repeats=n_repeats)
        if i % 10 == 0:
            df = pd.DataFrame.from_dict(results)
            df.to_csv("results_gpu.csv", index=False)

    df = pd.DataFrame.from_dict(results)

    df.to_csv("results_gpu.csv", index=False)
