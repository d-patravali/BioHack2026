from pathlib import Path
import os
from dotenv import load_dotenv
from openai import OpenAI


def main():
    # Load API key from .env
    load_dotenv(Path(__file__).resolve().parents[2] / ".env")

    client = OpenAI(
        api_key=os.environ["NEBIUS_API_KEY"],
        base_url="https://api.tokenfactory.nebius.com/v1/"
    )

    # Retrieve available models
    models = client.models.list()

    for model in models.data:
        print(model.id)


if __name__ == "__main__":
    main()
