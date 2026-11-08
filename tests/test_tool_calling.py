
import os
import json
from pathlib import Path
from dotenv import load_dotenv
from openai import OpenAI

# Load environment variables
PROJECT_ROOT = Path(__file__).resolve().parent.parent
load_dotenv(PROJECT_ROOT / ".env")

client = OpenAI(
    api_key=os.environ["NEBIUS_API_KEY"],
    base_url="https://api.tokenfactory.nebius.com/v1/"
)

MODEL = os.environ["NEBIUS_MODEL"]


# 1. Define a Python tool
def get_subtypes(disease: str):
    """Return supported disease subtypes."""

    subtypes = {
        "colorectal_cancer": [
            "GI.CIN",
            "GI.HM-indel",
            "GI.HM-SNV",
            "GI.GS"
        ]
    }

    return subtypes.get(disease, [])


# 2. Describe the tool to the LLM
tools = [
    {
        "type": "function",
        "function": {
            "name": "get_subtypes",
            "description": "Retrieve available subtypes for a disease.",
            "parameters": {
                "type": "object",
                "properties": {
                    "disease": {
                        "type": "string",
                        "description": "Disease identifier, e.g. colorectal_cancer"
                    }
                },
                "required": ["disease"]
            }
        }
    }
]


# 3. Ask the model a question
messages = [
    {
        "role": "system",
        "content": (
            "You are a biomedical research agent. "
            "Use available tools to retrieve factual data. "
            "Never invent supported diseases or subtypes."
        )
    },
    {
        "role": "user",
        "content": "What colorectal cancer subtypes are available?"
    }
]

response = client.chat.completions.create(
    model=MODEL,
    messages=messages,
    tools=tools,
    tool_choice="auto"
)

assistant_message = response.choices[0].message

# 4. Check whether the LLM requested a tool
if assistant_message.tool_calls:

    messages.append(assistant_message)

    for tool_call in assistant_message.tool_calls:

        print("Tool requested:", tool_call.function.name)
        print("Arguments:", tool_call.function.arguments)

        arguments = json.loads(tool_call.function.arguments)

        if tool_call.function.name == "get_subtypes":
            result = get_subtypes(**arguments)
        else:
            raise ValueError("Unknown tool requested")

        print("Tool result:", result)

        # Return the tool result to the model
        messages.append({
            "role": "tool",
            "tool_call_id": tool_call.id,
            "content": json.dumps(result)
        })

    # 5. Let the model generate its final answer
    final_response = client.chat.completions.create(
        model=MODEL,
        messages=messages
    )

    print("\nFinal answer:")
    print(final_response.choices[0].message.content)

else:
    print("The model did not request a tool.")
    print(assistant_message.content)
