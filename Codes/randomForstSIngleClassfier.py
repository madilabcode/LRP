import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import scipy.stats as sc
import pickle


def test_sig_as_classfier(exp, up, down=list()):
    cls = np.intersect1d(up + down + ["ident"], exp.columns)
    sub_exp = exp[cls]

    up_ls, down_ls = random_forst_expresstion(sub_exp, "sig_compere")
    pup = len(np.intersect1d(up, up_ls["feture"])) / len(up_ls)
    pdown = len(np.intersect1d(down, down_ls["feture"])) / len(down_ls)
    print("pup = " + str(pup))
    print("pdown = " + str(pdown))
    return up_ls, down_ls


def random_forst_expresstion(exp, plot_name="result", returnModle=False):
    new_exp = exp.sample(frac=1)
    new_exp.index = [i for i in range(len(new_exp))]
    pred_list = list()
    prob_list = list()
    target_list = list()
    size = math.ceil(len(new_exp) / 10)
    clf_save = dict()
    feture_import = dict()
    imp_sort = list()
    imp1_pval = list()
    imp1_feture = list()
    imp1_stats = list()
    imp2_pval = list()
    imp2_feture = list()
    imp2_stats = list()
    imp1_imp = list()
    imp2_imp = list()

    for i in range(10):
        test_sub = new_exp[i * size:min((i + 1) * size, len(new_exp))]
        traning_sub = new_exp.drop(axis=0, labels=[j for j in range(i * size, min((i + 1) * size, len(new_exp)))])
        clf = RandomForestClassifier(n_estimators=500, max_features="sqrt", oob_score=True)
        # print(traning_sub)
        # trining_sub = sub_df.drop(axis = 1, lbels = "target")
        clf.fit(traning_sub.drop(axis=1, labels="ident"), traning_sub["ident"])
        clf_save[i] = clf
        pred_list.append(clf.predict(test_sub.drop(axis=1, labels="ident")))
        prob_list.append(clf.predict_proba(test_sub.drop(axis=1, labels="ident")))
        target_list.append(test_sub["ident"])

    index = print_output_of_random_forst(target_list, pred_list)
    roc_test(prob_list, target_list, plot_name)

    clf = clf_save[index]
    imp = clf.feature_importances_

    for index, col in enumerate(new_exp.drop(axis=1, labels="ident").columns):
        feture_import[imp[index]] = col

    for key in feture_import.keys():
        imp_sort.append(key)

    imp_sort.sort(reverse=True)

    for key in imp_sort:
        feture = feture_import[key]
        s1, p1 = ttest_feture(feture, new_exp)
        if s1 > 0 and p1 <= 0.05:
            imp1_pval.append(p1)
            imp1_feture.append(feture)
            imp1_imp.append(key)
            imp1_stats.append(s1)

        elif s1 < 0 and p1 <= 0.05:
            imp2_pval.append(p1)
            imp2_feture.append(feture)
            imp2_imp.append(key)
            imp2_stats.append(s1)

    dict1 = {"feture": imp1_feture, "importances value": imp1_imp, "P_value": imp1_pval, "statstic": imp1_stats}
    imp_feture_up = pd.DataFrame(dict1)
    dict2 = {"feture": imp2_feture, "importances value": imp2_imp, "P_value": imp2_pval, "statstic": imp2_stats}
    imp_feture_down = pd.DataFrame(dict2)

    if returnModle:
        return imp_feture_up, imp_feture_down, clf

    return imp_feture_up, imp_feture_down


def use_ex_modle(mod, test_data, plot_name="result"):
    test_data = test_data.sample(frac=1)
    test_data.index = [i for i in range(len(test_data))]
    feture_import = dict()
    imp_sort = list()
    imp1_pval = list()
    imp1_feture = list()
    imp1_stats = list()
    imp2_pval = list()
    imp2_feture = list()
    imp2_stats = list()
    imp1_imp = list()
    imp2_imp = list()

    predict = mod.predict(test_data.drop(axis=1, labels="ident"))
    predict_prob = mod.predict_proba(test_data.drop(axis=1, labels="ident"))
    auc = roc_auc_score(test_data["ident"], predict_prob[:, 1])
    print("total auc: " + str(auc))
    plot_roc_carve(test_data["ident"], predict_prob[:, 1], auc, plot_name)
    imp = mod.feature_importances_

    for index, col in enumerate(test_data.drop(axis=1, labels="ident").columns):
        feture_import[imp[index]] = col

    for key in feture_import.keys():
        imp_sort.append(key)

    imp_sort.sort(reverse=True)

    for key in imp_sort:
        feture = feture_import[key]
        s1, p1 = ttest_feture(feture, test_data)
        if s1 > 0 and p1 <= 0.05:
            imp1_pval.append(p1)
            imp1_feture.append(feture)
            imp1_imp.append(key)
            imp1_stats.append(s1)

        elif s1 < 0 and p1 <= 0.05:
            imp2_pval.append(p1)
            imp2_feture.append(feture)
            imp2_imp.append(key)
            imp2_stats.append(s1)

    dict1 = {"feture": imp1_feture, "importances value": imp1_imp, "P_value": imp1_pval, "statstic": imp1_stats}
    imp_feture_up = pd.DataFrame(dict1)
    dict2 = {"feture": imp2_feture, "importances value": imp2_imp, "P_value": imp2_pval, "statstic": imp2_stats}
    imp_feture_down = pd.DataFrame(dict2)

    return imp_feture_up, imp_feture_down


def ttest_feture(feture, data):
    pop1 = data.loc[data["ident"] == 0, feture]
    pop2 = data.loc[data["ident"] == 1, feture]

    s1, p1 = sc.ttest_ind(np.array(pop1), np.array(pop2))

    return s1, p1


def print_output_of_random_forst(target_list, pred_list):
    n = len(target_list)
    total_correct = 0
    n_of_elements = 0
    max_p = 0
    max_index = 0
    for i in range(n):
        unite_arry = np.array(target_list[i]) + np.array(pred_list[i])
        unite_arry = unite_arry[unite_arry != 1]
        sub_correct = len(unite_arry)
        total_correct += sub_correct
        n_of_elements += len(target_list[i])

        if (total_correct / n_of_elements) > max_p:
            max_p = (total_correct / n_of_elements)
            max_index = i

        print("the " + str(i) + "test sub")
        print(np.array(target_list[i]))
        print(pred_list[i])
        print("in this sub test the correct % is  " + str(sub_correct / len(target_list[i])))
    print("total correct % is  " + str(total_correct / n_of_elements) + "\n")

    return max_index


def roc_test(prob_list, target_list, plot_name="result"):
    all_prob = []
    all_target = []

    for i in range(len(prob_list)):
        all_prob = np.concatenate((all_prob, np.array(prob_list[i][:, 1])))
        all_target = np.concatenate((all_target, np.array(target_list[i])))

    auc = roc_auc_score(all_target, all_prob)
    print("total auc: " + str(auc))

    fpr, tpr, threshold = roc_curve(all_target, all_prob)
    plt.plot(fpr, tpr, linestyle='--', label='Random Forest AUC')
    plt.plot([0, 1], [0, 1], "--")
    plt.title("Random Forest Classifier AUC Value: " + str(round(auc, 3)))
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.xlabel("FPR RATE")
    plt.ylabel("TPR RATE")
    plt.savefig(plot_name + ".pdf", bbox_inches='tight', transparent=True, pad_inches=0.5)

    return tpr, fpr


def plot_roc_carve(target, prob, auc, plot_Name="result"):
    fpr, tpr, threshold = roc_curve(target, prob)
    plt.plot(fpr, tpr, linestyle='--', label='Random Forest AUC')
    plt.plot([0, 1], [0, 1], "--")
    plt.title("Random Forest Classifier AUC Value: " + str(round(auc, 3)))
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.xlabel("FPR RATE")
    plt.ylabel("TPR RATE")
    plt.savefig(plot_Name + ".pdf", bbox_inches='tight', transparent=True, pad_inches=0.5)
    return fpr, tpr


if __name__ == "__main__":
    pass
    # exp = pd.read_csv("exp.csv", index_col=None)
# print(exp)
# l1,l2 = random_forst_expresstion(exp)
