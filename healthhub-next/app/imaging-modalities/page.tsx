import React from "react";

import ServiceCard from "../../components/services/ServiceCard";
import {Divider} from "@nextui-org/divider";
/**
 * An array of service objects, each representing a medical imaging modality.
 * Each service object contains the following properties:
 *
 * @property {string} title - The title of the imaging modality.
 * @property {string} imagePath - The path to the thumbnail image representing the modality.
 * @property {string} description - A brief description of the imaging modality and its AI capabilities.
 * @property {string} path - The URL path to the detailed page of the imaging modality.
 */
const services = [
  {
    title: "Dibatic Retinopathy Detection",
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
  {
    title: "COVID-19 Detection",
    imagePath: "/modalities_thumbnails/covid19.jpeg",
    description:
      "AI-assisted detection of COVID-19 disease.",
    path: "/imaging-modalities/covid-19-prediction",
  },
  {
    title: "Image Processing",
    imagePath: "/modalities_thumbnails/image-proccesing.jpeg",
    description:
      "Standard Image Proceesing Tool",
    path: "/imaging-modalities/image-processing",
  },
  {
    title: "Invasive Ductal Carcinoma Detection",
    imagePath: "/modalities_thumbnails/breast-cancer.jpg",
    description:
      "AI-powered detection of invasive ductal carcinoma in breast cancer.",
    path: "/imaging-modalities/invasive-ductal-carcinoma",
  },
  {
    title: "Skin Cancer Detection",
    imagePath: "/modalities_thumbnails/breast-cancer.jpg",
    description:
      "Detect Melignant and Benign cancer in dermascopy images.",
    path: "/imaging-modalities/skin-cancer",
  }
];
// const services2 = [{
//   title: "Image Processing",
//   imagePath: "/modalities_thumbnails/image-proccesing.jpeg",
//   description:
//     "Standard Image Proceesing Tool",
//   path: "/imaging-modalities/image-processing",
// },]
/**
 * Renders the Services component which displays a list of AI-powered services.
 *
 * @returns {JSX.Element} The rendered Services component.
 *
 * The component includes:
 * - A main container with a maximum width and padding.
 * - A heading with the title "Our AI-Powered Services".
 * - A paragraph with a brief description.
 * - A grid layout that displays a list of service cards.
 *
 * Each service card is rendered using the `ServiceCard` component and includes:
 * - `title`: The title of the service.
 * - `imagePath`: The path to the service's image.
 * - `description`: A brief description of the service.
 * - `path`: The path to the service's detailed page.
 */
export default function Services() {
  return (
    <div className="max-w-7xl mx-auto p-6">
      <h1 className="text-4xl font-extrabold text-center mb-8">
        Our AI-Powered Services
      </h1>
      <p className="text-lg text-gray-600 text-center mb-12">
        Explore our advanced tools and services. Click on a service to learn
        more.
      </p>

      <div className="grid gap-6 grid-cols-1 sm:grid-cols-2 md:grid-cols-3">
        {services.map((service, index) => (
          <ServiceCard
            key={index}
            description={service.description}
            imagePath={service.imagePath}
            path={service.path}
            title={service.title}
          />
        ))}
      </div>
      <Divider className='my-4'/>
      <div className="text-center mt-8">
        <h2 className="text-2xl font-bold mb-4">Image Proccesing</h2>
        <p className="text-lg text-gray-600">

        </p>
        <div className="grid gap-6 grid-cols-1 sm:grid-cols-2 md:grid-cols-3">
        {/* {services2.map((service, index) => (
          <ServiceCard
            key={index}
            description={service.description}
            imagePath={service.imagePath}
            path={service.path}
            title={service.title}
          />
        ))} */}
      </div>
        </div>

    </div>
  );
}
