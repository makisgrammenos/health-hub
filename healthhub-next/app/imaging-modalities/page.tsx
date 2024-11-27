import React from "react";
import ServiceCard from "../../components/services/ServiceCard";

const services = [
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

export default function Services() {
  return (
    <div className="max-w-7xl mx-auto p-6">
      <h1 className="text-4xl font-extrabold text-center mb-8">
        Our AI-Powered Services
      </h1>
      <p className="text-lg text-gray-600 text-center mb-12">
        Explore our advanced tools and services. Click on a service to learn more.
      </p>

      <div className="grid gap-6 grid-cols-1 sm:grid-cols-2 md:grid-cols-3">
        {services.map((service, index) => (
          <ServiceCard
            key={index}
            title={service.title}
            imagePath={service.imagePath}
            description={service.description}
            path={service.path}
          />
        ))}
      </div>
    </div>
  );
}
