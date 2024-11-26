import React from "react";
import NextLink from "next/link";

const modalities = [
  {
    title: "Fundus Imaging",
    imagePath: "/modalities_thumbnails/diabetic_retinopathy.jpg",
    description: "Leverage AI to detect diabetic retinopathy with precision.",
    path: "/imaging-modalities/diabetic-retinopathy",
  },
  {
    title: "Brain Tumor Segmentation",
    imagePath: "/modalities_thumbnails/brain_seg.jpg",
    description:
      "AI-powered 3D MRI segmentation for accurate brain tumor analysis.",
    path: "/imaging-modalities/brain-tumor-segmentation",
  },
  {
    title: "Chest Pathology Classification",
    imagePath: "/modalities_thumbnails/chest.jpg",
    description:
      "AI-assisted detection of chest pathologies using X-ray imaging.",
    path: "/imaging-modalities/chest-pathology-classification",
  },
];

export default function ImagingModalities() {
  return (
    <div className="max-w-7xl mx-auto p-6">
      <h1 className="text-4xl font-extrabold text-center mb-8">
        AI-Powered Image Analysis Tools
      </h1>
      <p className="text-lg text-gray-600 text-center mb-12">
        Discover advanced medical imaging tools and analysis. Click on a modality to learn more.
      </p>

      <div className="grid gap-6 grid-cols-1 sm:grid-cols-2 md:grid-cols-3">
        {modalities.map((modality, index) => (
          <div
            key={index}
            className="group relative cursor-pointer rounded-lg overflow-hidden shadow-lg hover:shadow-2xl transition-shadow"
          >
            <NextLink href={modality.path}>
              <div className="relative w-full h-60">
                <img
                  src={modality.imagePath}
                  alt={modality.title}
                  className="object-cover w-full h-full"
                />
              </div>
              <div className="absolute inset-0 bg-black bg-opacity-50 group-hover:bg-opacity-30 transition-all"></div>
              <div className="absolute bottom-4 left-4">
                <h2 className="text-lg font-bold text-white">{modality.title}</h2>
                <p className="text-sm text-gray-300">{modality.description}</p>
              </div>
            </NextLink>
          </div>
        ))}
      </div>
    </div>
  );
}
