from xgboost.sklearn import XGBClassifier
import joblib
from trainIO import *
from sklearn.metrics import f1_score, accuracy_score, recall_score, precision_score, roc_auc_score
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import RFE


def custom_scorer(estimator, y_true, y_pred):
    return f1_score(y_true, y_pred)

if __name__ == "__main__":
    print('XGBM Start: ' + str(datetime.datetime.now()))
    f_names = []

    with open('./learn.window.selection', 'r') as f:
        for line in f:
            f_names.append(line.strip())
    f_names.pop(0)
    train_x, train_y, train_name, sc = loadfile('./learn.vector.train')
    xgb_model = XGBClassifier()
    # define the grid search parameters
    param_grid = {
        'learning_rate': [0.05],
        'n_estimators': [170],
        'max_depth': [3],
        'min_child_weight': [1],
        'gamma': [0.4],
        'colsample_bytree': [0.9],
        'reg_alpha': [0.5],
        'reg_lambda': [0.8],
        'n_jobs': [4],
        'nthread': [4],
    }

    grid_search = GridSearchCV(estimator=xgb_model,
                               param_grid=param_grid,
                               scoring=['accuracy', 'roc_auc'],
                               cv=5,
                               verbose=0,
                               n_jobs=5,
                               refit="roc_auc"
                               )
    grid_search.fit(train_x, train_y)

    print("Best parameters: ", grid_search.best_params_)
    print("Best score: ", grid_search.best_score_)

    xgbm = grid_search.best_estimator_
    train_x_new = train_x
    total = train_x.shape[1]
    max_score = 0
    max_n_f = 0
    n_select = total - 1
    removed_features = []
    while n_select > 2:
        step = max(math.ceil(n_select * 0.03), 2)
        selector = RFE(xgbm, n_features_to_select=n_select - step, verbose=0, step=step)
        train_x_new = selector.fit_transform(train_x_new, train_y)
        grid_search.fit(train_x_new, train_y)
        f_names_new = [value for value, flag in zip(f_names, selector.support_) if flag]
        # get removed feature by diff f_names and f_names_new, then append to removed_features
        removed_feature = [value for value in f_names if value not in f_names_new]
        removed_features.extend(removed_feature)
        f_names = f_names_new

        xgbm = grid_search.best_estimator_
        if (max_score <= grid_search.best_score_):
            max_score = grid_search.best_score_
            f_list = f_names
            max_n_f = n_select - step
            model = xgbm
            print(f"new Max Score: {max_score} at {max_n_f} features.")
        print(f"Score: {grid_search.best_score_} at {n_select - step} features.")
        n_select = n_select - step

    print(f"Max Score: {max_score} at {max_n_f} features.")
    print(f"feature:", ','.join(f_list))
    lines_with_newlines = [line + '\n' for line in f_list]
    removed_features.extend(f_names)
    with open('output.RFE.txt', 'w', encoding='utf-8') as file:
        file.writelines(lines_with_newlines)
    removed_features_with_newlines = [line + '\n' for line in removed_features]
    with open('output.RFE.full.txt', 'w', encoding='utf-8') as file:
        file.writelines(removed_features_with_newlines)
    model.save_model('./xgbm.model.json')
    joblib.dump(sc, './xgbm.sc.pkl')

    print('End: ' + str(datetime.datetime.now()))
