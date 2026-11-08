import os
from dotenv import load_dotenv
from openai import OpenAI

load_dotenv()

client = OpenAI(
    api_key=os.environ["NEBIUS_API_KEY"],
    base_url="https://api.tokenfactory.nebius.com/v1/"
)

response = client.chat.completions.create(
    model=os.environ["NEBIUS_MODEL"],
    messages=[
        {
            "role": "user",
            "content": "What is a graph neural network? Answer in one sentence."
        }
    ]
)

print(response.choices[0].message.content)