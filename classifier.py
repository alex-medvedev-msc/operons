import pandas
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn import svm
from sklearn.metrics import log_loss


def load_features(path):
    values = pandas.read_csv(path).values
    return train_test_split(values[:, :-1], values[:, -1], test_size=0.25)


def main():
    path = "features.csv"

    X_train, X_test, y_train, y_test = load_features(path)
    gbt_classifier = GradientBoostingClassifier()
    gbt_classifier.fit(X_train, y_train)
    y_pred = gbt_classifier.predict_proba(X_test)[:, 1]
    preds = pandas.DataFrame.from_dict({"pred": list(y_pred), "test": list(y_test)})
    print(preds)
    print(log_loss(y_test, y_pred))

    svm_classifier = svm.SVC(probability=True)
    svm_classifier.fit(X_train, y_train)

    y_pred = svm_classifier.predict_proba(X_test)[:, 1]
    preds = pandas.DataFrame.from_dict({"pred": list(y_pred), "test": list(y_test)})
    print(preds)
    print(log_loss(y_test, y_pred))



if __name__ == '__main__':
    main()