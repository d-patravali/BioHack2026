import os
from dotenv import load_dotenv
from openai import OpenAI

# Load API key from .env
load_dotenv()

client = OpenAI(
    api_key=os.environ["NEBIUS_API_KEY"],
    base_url="https://api.tokenfactory.nebius.com/v1/"
)

# Retrieve available models
models = client.models.list()

for model in models.data:
    print(model.id)