from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse, StreamingResponse
import pandas as pd
from catboost import CatBoostClassifier
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
import matplotlib.pyplot as plt
import seaborn as sns
import io

router = APIRouter()

def rank_list(lst):
    return {elm: len(lst) - 1 - i for i, elm in enumerate(lst)}

def borda_aggregation(loflists):
    list_ranks = [rank_list(l) for l in loflists]
    feature_set = {i for i in [el for nl in loflists for el in nl]}
    return {e: sum([lr.get(e, 0) for lr in list_ranks]) for e in feature_set}

def create_sorted_df(result):
    df = pd.DataFrame(list(result.items()), columns=["Feature", "Borda Rank"])
    df.sort_values(by="Borda Rank", ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df

def borda_df(loflists):
    borda_dict = borda_aggregation(loflists)
    return create_sorted_df(borda_dict)

def feature_importance_df_creation(model, model_label, feature_names):
    importance_df = pd.DataFrame(model.feature_importances_)
    importance_df = importance_df.rename(columns={0: model_label})
    importance_df["feature_name"] = feature_names
    importance_df = importance_df.sort_values(by=model_label, ascending=False)
    return importance_df

@router.get("/feature_selection")
def feature_selection():
    try:
        # Load the dataset
        df = pd.read_csv("routers/clinical/uci_hf_df.csv")
        tag = pd.DataFrame(df, columns=["tag"])
        df = df.drop(columns=['tag'])

        # Initialize models
        xgbclass = XGBClassifier(eval_metric='logloss', verbosity=0, importance_type='gain', silent=True)
        cbc = CatBoostClassifier(iterations=100, verbose=0)
        lgbm = LGBMClassifier(importance_type='split')

        # Train models
        xgbclass.fit(df, tag["tag"])
        cbc.fit(df, tag["tag"])
        lgbm.fit(df, tag["tag"])

        # Feature importance
        xgb_importance = feature_importance_df_creation(xgbclass, "XgBoost_importance", list(df.columns))
        cbc_importance = feature_importance_df_creation(cbc, "CatBoost_importance", list(df.columns))
        lgbm_importance = feature_importance_df_creation(lgbm, "LightGBM_importance", list(df.columns))

        # Prepare lists for Borda rank
        xgboost_imp_lst = list(xgb_importance["feature_name"])
        catboost_imp_lst = list(cbc_importance["feature_name"])
        lgbm_imp_lst = list(lgbm_importance["feature_name"])

        # Compute Borda results
        borda_results = borda_df([xgboost_imp_lst, catboost_imp_lst, lgbm_imp_lst])

        # Create the plot
        plt.figure(figsize=(8, 6))
        ax = sns.scatterplot(data=borda_results, x="Borda Rank", y="Feature", s=100, color="steelblue", label=" ")
        x_start = borda_results["Borda Rank"].min()
        ax.set_xlim(left=x_start)
        ax.set_facecolor("whitesmoke")
        for _, row in borda_results.iterrows():
            ax.plot([x_start, row["Borda Rank"]], [row["Feature"], row["Feature"]], color="steelblue", linewidth=0.5)
        ax.set_title("Borda Consensus Feature Importance")
        ax.set_xlabel("Importance")
        ax.set_ylabel("Feature")
        ax.legend_.remove()
        plt.tight_layout()

        # Save the plot to a buffer
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        plt.close()
        buf.seek(0)

        return StreamingResponse(buf, media_type="image/png")

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error in feature selection: {str(e)}")

@router.get("/borda_results")
def get_borda_results():
    try:
        # Load the dataset
        df = pd.read_csv("uci_hf_df.csv")
        tag = pd.DataFrame(df, columns=["tag"])
        df = df.drop(columns=['tag'])

        # Initialize models
        xgbclass = XGBClassifier(eval_metric='logloss', verbosity=0, importance_type='gain', silent=True)
        cbc = CatBoostClassifier(iterations=100, verbose=0)
        lgbm = LGBMClassifier(importance_type='split')

        # Train models
        xgbclass.fit(df, tag["tag"])
        cbc.fit(df, tag["tag"])
        lgbm.fit(df, tag["tag"])

        # Feature importance
        xgb_importance = feature_importance_df_creation(xgbclass, "XgBoost_importance", list(df.columns))
        cbc_importance = feature_importance_df_creation(cbc, "CatBoost_importance", list(df.columns))
        lgbm_importance = feature_importance_df_creation(lgbm, "LightGBM_importance", list(df.columns))

        # Prepare lists for Borda rank
        xgboost_imp_lst = list(xgb_importance["feature_name"])
        catboost_imp_lst = list(cbc_importance["feature_name"])
        lgbm_imp_lst = list(lgbm_importance["feature_name"])

        # Compute Borda results
        borda_results = borda_df([xgboost_imp_lst, catboost_imp_lst, lgbm_imp_lst])

        # Return the results as JSON
        return JSONResponse(content=borda_results.to_dict(orient="records"))

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error computing Borda results: {str(e)}")
