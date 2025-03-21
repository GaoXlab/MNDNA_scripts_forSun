import joblib
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import LabelEncoder
from xgboost.sklearn import XGBClassifier

from trainIO import *


def top_k_accuracy_per_label(y_true, y_pred_proba, k=2):
    # Get the unique labels
    labels = np.unique(y_true)

    # Initialize an empty dictionary to store the top-k accuracy for each label
    label_acc_dict = {}

    for label in labels:
        # Filter out the specific label
        indices = (y_true == label)
        y_true_label = y_true[indices]
        y_pred_proba_label = y_pred_proba[indices]

        # Get the top k predictions
        top_k_pred = np.argsort(y_pred_proba_label, axis=1)[:, -k:]

        # Check if the true labels are in top k predictions
        match_array = np.any(top_k_pred == y_true_label[:, None], axis=1)

        # Calculate top k accuracy
        top_k_acc = np.mean(match_array)

        # Add the top-k accuracy to the dictionary
        label_acc_dict[label] = top_k_acc

    return label_acc_dict


def print_cm(cm, le):
    class_labels = le.inverse_transform([i for i in range(len(le.classes_))])
    print(',' + ','.join(str(label) for label in class_labels))
    for encoded_label, row in enumerate(cm):
        original_label = le.inverse_transform([encoded_label])[0]
        # calculate（True Positives）
        tp = cm[encoded_label, encoded_label]
        # calculate（Predicted Positives）
        pp = np.sum(cm[encoded_label, :])
        # calculate（Accuracy）
        accuracy = tp / pp if pp > 0 else 0
        print(f"{original_label}," + ','.join(str(value) for value in row) + f",{accuracy:.2f}")


if __name__ == "__main__":
    print('XGBM Start: ' + str(datetime.datetime.now()))

    train_x, train_y, train_name, sc = loadfile('./learn.vector.train')
    test_x, test_y, test_name, sc = loadfile('./learn.vector.test', sc)
    le = LabelEncoder()
    train_y = le.fit_transform(train_y)
    test_y = le.transform(test_y)

    # xgboost model
    xgbm = XGBClassifier(
        objective='multi:softprob',
        num_class=len(set(train_y)),
        learning_rate=0.1, n_estimators=350, max_depth=3, min_child_weight=1,
        seed=1234, subsample=0.8, colsample_bytree=1, gamma=0.7, reg_alpha=1, reg_lambda=1, n_jobs=12)

    removed_feature = []
    fit2 = xgbm.fit(train_x, train_y)
    xgbm_predict = fit2.predict_proba(test_x)

    fit2.save_model('./xgbm.model.json')
    joblib.dump(sc, './xgbm.sc.pkl')

    xgbm_predict_train = fit2.predict(train_x)
    xgbm_predict_train_label = le.inverse_transform(xgbm_predict_train)
    xgbm_predict_prob_train = fit2.predict_proba(train_x)

    y_pred = fit2.predict(test_x)

    cm = confusion_matrix(test_y, y_pred)
    train_cm = confusion_matrix(train_y, xgbm_predict_train)

    print("Confusion Matrix for training set:")
    print_cm(train_cm, le)
    print("Confusion Matrix for test set:")
    print_cm(cm, le)
    print("TOP-2 accuracy per label:")
    print(top_k_accuracy_per_label(test_y, xgbm_predict, k=2))

    # shuffle test label and calculate cm again
    print("Confusion Matrix for test set after shuffling:")
    shuffle_test_y = test_y.copy()
    np.random.seed(1234)
    np.random.shuffle(shuffle_test_y,)
    shuffled_cm = confusion_matrix(shuffle_test_y, y_pred)
    print_cm(shuffled_cm, le)
    print(top_k_accuracy_per_label(shuffle_test_y, xgbm_predict, k=2))

    # output each test sample's pred result to file './xgbm.predict_result'
    y_pred_label = le.inverse_transform(y_pred)
    y_test_label = le.inverse_transform(test_y)
    with open('./xgbm.predict_result', 'w') as file:
        for i in range(len(test_y)):
            file.write(f"{test_name[i]},{y_pred_label[i]},{y_test_label[i]}\n")
    # output each test sample's pred prob to file './xgbm.predict_prob_result'
    with open('./xgbm.predict_prob_result', 'w') as file:
        for i in range(len(xgbm_predict)):
            file.write(f"{test_name[i]},{','.join(str(value) for value in xgbm_predict[i])},{y_pred_label[i]}\n")
    # output each train sample's pred prob to file './xgbm.predict_prob_train'
    with open('./xgbm.predict_prob_train', 'w') as file:
        for i in range(len(xgbm_predict_prob_train)):
            file.write(f"{train_name[i]},"
                       + f"{','.join(str(value) for value in xgbm_predict_prob_train[i])},{xgbm_predict_train_label[i]}\n")

    #save shuffled test label and pred prob to file
    with open('./xgbm.predict_prob_shuffled', 'w') as file:
        for i in range(len(shuffle_test_y)):
            file.write(f"{test_name[i]},{','.join(str(value) for value in xgbm_predict[i])},{le.inverse_transform([shuffle_test_y[i]])[0]}\n")