'use client';

import React from "react";
import { Card, CardBody,CardFooter,CardHeader} from "@nextui-org/card";
import {Divider} from "@nextui-org/divider";
import { Image } from "@nextui-org/image";
const modalities = [
  {
    title: "Fundus Imaging",
    imagePath: "/modalities_thumbnails/diabetic_retinopathy.jpg",
    description: "Diabetic Retinopathy Detection"
  },
  {
    title: "Brain Tumor Segmentation",
    imagePath: "/modalities_thumbnails/brain_seg.jpg",
    description: "3D MRI Brain Tumor Segmentation"
  },
  {
    title: "Chest Pathology Classification",
    imagePath: "/modalities_thumbnails/chest.jpg",
    description: "Chest Pathology Classification"
  }
];

export default function ImagingModalities() {
  return (
    <div className="gap-6 grid grid-cols-1 sm:grid-cols-3 max-w-7xl mx-auto p-6">
      {modalities.map((modality, index) => (
        <Card
          key={index}
          shadow="sm"
          isPressable
          onPress={() => console.log(`Navigating to ${modality.title}`)}
        >
          <CardHeader className="text-center">{modality.title}</CardHeader>
          <Divider/>
          <CardBody className="overflow-visible p-0">

            
            <Image
              shadow="sm"
              radius="lg"
              width="100%"
              alt={modality.title}
              className="w-full object-cover h-[200px]"
              src={modality.imagePath}
            />
          </CardBody>
          <CardFooter className="text-small justify-between">
            {/* <b>{modality.title}</b> */}
            <p className="text-default-500">{modality.description}</p>
          </CardFooter>
        </Card>
      ))}
    </div>
  );
}
