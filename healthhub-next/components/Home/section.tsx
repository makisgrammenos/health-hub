import { Button } from "@nextui-org/button";
import Image from "next/image";

type SectionProps = {
  title: string;
  description: string;
  buttonText: string;
  buttonLink: string;
  imageSrc: string;
};

const Section: React.FC<SectionProps> = ({
  title,
  description,
  buttonText,
  buttonLink,
  imageSrc,
}) => {
  return (
    <div className="relative w-screen h-[400px] group">
      <Image
        src={imageSrc}
        alt={`${title} Background`}
        layout="fill"
        objectFit="cover"
        className="z-0 transition-transform duration-500 ease-in-out group-hover:scale-105 group-hover:blur-sm"
      />
      <div className="absolute inset-0 bg-gradient-to-t from-black/70 via-black/40 to-transparent" />
      <div className="absolute left-8 top-1/2 transform -translate-y-1/2 text-left text-white px-4">
        <h2 className="text-4xl font-bold">{title}</h2>
        <p className="mt-4 text-lg">{description}</p>
        <Button
          as="a"
          href={buttonLink}
          color="primary"
          radius="lg"
          size="md"
          className="mt-6"
        >
          {buttonText}
        </Button>
      </div>
    </div>
  );
};

export default Section;
