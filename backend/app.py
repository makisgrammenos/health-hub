from fastapi import FastAPI
from routers import diabetic_retinopathy,brain_tumor_segmentation,chestXray,image_processing, covid_prediction,breast_cancer,skin_cancer,region_of_interset
from fastapi.middleware.cors import CORSMiddleware
from routers.clinical import preprocessing

# Define allowed origins
origins = [
    "http://localhost:3000",  # Next.js development server
    # "https://your-nextjs-domain.com",  # Production domain
]
app = FastAPI(
    title="Medical Imaging Prediction API",
    description="A scalable API for multiple medical imaging models.",
    version="1.0.0",
)
# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],  # Allow all HTTP methods
    allow_headers=["*"],  # Allow all headers
)


# Include routers for each model
app.include_router(diabetic_retinopathy.router, prefix="/imaging/diabetic-retinopathy", tags=["Diabetic Retinopathy"])
app.include_router(brain_tumor_segmentation.router, prefix="/imaging/brain-tumor", tags=["Brain Tumor Segmentation"])
app.include_router(chestXray.router, prefix="/imaging/chest-x-ray", tags=["Chest X-Ray"])
app.include_router(image_processing.router, prefix="/imaging/image-processing", tags=["Image Processing"])
app.include_router(covid_prediction.router, prefix="/imaging/covid", tags=["COVID-19 Prediction"])
app.include_router(breast_cancer.router, prefix="/imaging/breast-cancer", tags=["Breast Cancer Prediction"])
app.include_router(skin_cancer.router, prefix="/imaging/skin-cancer", tags=["Skin Cancer Prediction"])
app.include_router(region_of_interset.router, prefix="/imaging/roi", tags=["Region of Interest"])


app.include_router(preprocessing.router, prefix="/clinical", tags=["Clinical Data Preprocessing"])

@app.get("/")
async def root():
    return {"message": "Welcome to the Medical Imaging Prediction API"}
