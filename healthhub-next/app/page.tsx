import Section from "../components/Home/section";

export default function Home() {
  return (
    <div className="w-full overflow-x-hidden">
      <Section
        title="Imaging Tools"
        description="Access advanced imaging tools for data analysis and visualization."
        buttonText="Explore Imaging Tools"
        buttonLink="/imaging-modalities"
        imageSrc="/banner/medical-imaging.jpg" // Replace with your actual image path
      />
      <Section
        title="Clinical Tools"
        description="Leverage clinical data analysis tools to derive actionable insights."
        buttonText="Explore Clinical Tools"
        buttonLink="/clinical"
        imageSrc="/banner/clinical.jpg" // Replace with your actual image path
      />
      <Section
        title="Molecular Tools"
        description="Dive deep into molecular-level data with powerful analysis tools."
        buttonText="Explore Molecular Tools"
        buttonLink="/molecular"
        imageSrc="/banner/molecular.png" // Replace with your actual image path
      />
    </div>
  );
}
