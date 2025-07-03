# SPDX-FileCopyrightText: The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import csv
import json
import os
import re
from datetime import datetime
from pathlib import Path

import requests
import typer
from tqdm import tqdm
from tqdm.utils import CallbackIOWrapper

app = typer.Typer()

ZENODO_API_URL = "https://sandbox.zenodo.org/api/deposit/depositions"
VERSIONS_CSV = Path("data/versions.csv")


def get_access_token() -> str:
    """
    Retrieve the Zenodo API access token from the environment or prompt the user.

    Returns
    -------
    str
        The Zenodo API access token.
    """
    token = os.environ.get("ZENODO_ACCESS_TOKEN")
    if not token:
        token = typer.prompt(
            "ZENODO_ACCESS_TOKEN not found in environment. You can set it up using 'export ZENODO_ACCESS_TOKEN='<token>'. "
            "Please enter your Zenodo API token"
        )
    return token


def read_versions_csv() -> list[dict]:
    """
    Read the versions CSV file and return its contents as a list of dictionaries.

    Returns
    -------
    list of dict
        List of rows from the versions CSV file.
    """
    if not VERSIONS_CSV.exists():
        raise FileNotFoundError(
            f"Versions CSV file not found at {VERSIONS_CSV}. Please create and populate it first."
        )
    rows = []
    with open(VERSIONS_CSV, newline="") as f:
        for row in csv.DictReader(f):
            # Fix tags column to always be a list
            row["tags"] = (
                json.loads(row["tags"].replace("'", '"')) if row["tags"] else []
            )
            rows.append(row)
    return rows


def write_versions_csv(rows: list[dict]):
    """
    Write a list of dictionaries to the versions CSV file.

    Parameters
    ----------
    rows : list of dict
        The rows to write to the CSV file.
    """

    for row in rows:
        # Ensure tags are always a JSON string with single quotes
        row["tags"] = json.dumps(row["tags"]).replace('"', "'")

    with open(VERSIONS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["dataset", "source", "version", "tags", "url", "note"],
            quotechar='"',
            quoting=csv.QUOTE_ALL,
        )
        writer.writeheader()
        writer.writerows(rows)


def prompt_choice(
    options: list[str], prompt_text: str, other_text: str | None = None
) -> str:
    """
    Prompt the user to select an option from a list by number,
    with optional support for a custom value ("other").

    Parameters
    ----------
    options : list of str
        The available options.
    prompt_text : str
        The prompt to display.
    other_text : str or None, optional
        The label for the custom value option (default is None, disables "other").

    Returns
    -------
    str
        The selected or custom value.
    """
    display_options = options.copy()
    if other_text is not None:
        display_options.append(other_text)

    # List options with prepended numbers
    typer.echo(prompt_text)
    for idx, option in enumerate(display_options, 1):
        typer.echo(f"{idx}: {option}")

    # Prompt user for choice
    if len(display_options) == 1:
        # Only one option, use as default
        choice = typer.prompt("Select (1)", default="1")
    else:
        choice = typer.prompt(f"Select (1-{len(display_options)})")

    # Validate choice
    if choice.isdigit() and 1 <= int(choice) <= len(display_options):
        selected = display_options[int(choice) - 1]
        if other_text is not None and selected == other_text:
            return typer.prompt("Please specify the new value")
        return selected
    typer.echo("Invalid selection. Please enter a number from the list.")

    return prompt_choice(options, prompt_text, other_text)


def get_latest_versions(rows, source):
    """
    Get the latest versions for a given source.

    Parameters
    ----------
    rows : list of dict
        The rows from the versions CSV.
    source : str
        The source to filter by.

    Returns
    -------
    list of dict
        The filtered and sorted latest versions.
    """
    latest = [
        row
        for row in rows
        if row["source"] == source and "latest" in row["tags"] and row.get("url")
    ]
    return sorted(latest, key=lambda r: r["dataset"])


def get_primary_datasets(rows):
    """
    Get the list of primary datasets not present in the archive.

    Parameters
    ----------
    rows : list of dict
        The rows from the versions CSV.

    Returns
    -------
    list of str
        The primary datasets.
    """
    primary = {row["dataset"] for row in rows if row["source"] == "primary"}
    archive = {row["dataset"] for row in rows if row["source"] == "archive"}
    return sorted(list(primary - archive))


def get_dataset_versions(rows, dataset):
    """
    Get all versions for a given dataset.

    Parameters
    ----------
    rows : list of dict
        The rows from the versions CSV.
    dataset : str
        The dataset name.

    Returns
    -------
    list of str
        The versions for the dataset.
    """
    return [row["version"] for row in rows if row["dataset"] == dataset]


def get_new_archive_folders(dataset, known_versions):
    """
    Get new archive folders for a dataset that are not in known versions.

    Parameters
    ----------
    dataset : str
        The dataset name.
    known_versions : list of str
        The known versions.

    Returns
    -------
    list of str
        The new archive folder names.
    """
    archive_path = Path(f"data/{dataset}/archive")
    if not archive_path.exists():
        return []
    return sorted(
        [
            f.name
            for f in archive_path.iterdir()
            if f.is_dir() and f.name not in known_versions
        ]
    )


def get_potential_datasets():
    """
    Get all dataset folders that have an archive subfolder and are thus potential candidates for a new upload.

    Returns
    -------
    list of str
        The dataset names that have an archive subfolder.
    """
    # Find all dataset folders in the 'data' directory;
    # indicated by the presence of an 'archive' subfolder with contents
    return sorted(set(f.parts[-3] for f in Path().rglob("data/*/archive/*/")))


def get_archive_folders(dataset):
    """
    Get all archive folders for a dataset.

    Parameters
    ----------
    dataset : str
        The dataset name.

    Returns
    -------
    list of str
        The archive folder names.
    """
    archive_path = Path(f"data/{dataset}/archive")
    if not archive_path.exists():
        return []
    return sorted([f.name for f in archive_path.iterdir() if f.is_dir()])


def create_zenodo_deposition(
    token: str, api_url: str, metadata: dict, files: list[Path]
) -> requests.Response:
    """
    Create a new Zenodo deposition with the specified metadata and files.

    Parameters
    ----------
    token : str
        The Zenodo API access token to use.
    api_url : str, optional
        The Zenodo API URL for depositions.
    metadata : dict, optional
        The metadata for the deposition.
    files : list of Path, optional
        The files to upload to the deposition.

    Returns
    -------
    Response
        The response from the Zenodo API after creating the deposition.
    """
    r = requests.post(
        api_url,
        params={"access_token": token},
        json={},
        headers={"Content-Type": "application/json"},
    )

    r.raise_for_status()
    bucket_url = r.json()["links"]["bucket"]
    deposition_url = r.json()["links"]["self"]

    # Update metadata of the deposition
    r = requests.put(
        deposition_url,
        params={"access_token": token},
        json={"metadata": metadata},
        headers={"Content-Type": "application/json"},
    )

    r.raise_for_status()

    # Upload files one-by-one to the bucket
    for file in files:
        with open(file, "rb") as fp:
            # Progressbar on request.puts https://gist.github.com/tyhoff/b757e6af83c1fd2b7b83057adf02c139
            with tqdm(
                total=file.stat().st_size,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=f"Uploading {file.name}",
            ) as pbar:
                wrapped_file = CallbackIOWrapper(pbar.update, fp, "read")
                upload_response = requests.put(
                    f"{bucket_url}/{file.name}",
                    data=wrapped_file,
                    params={"access_token": token},
                )
                upload_response.raise_for_status()
                # Add done to the progress bar
                pbar.update(file.stat().st_size)

    return r


def publish_zenodo_deposition(token: str, deposition_id: int) -> requests.Response:
    """
    Publish a Zenodo deposition.

    Parameters
    ----------
    token : str
        The Zenodo API access token.
    deposition_id : int
        The ID of the deposition to publish.

    Returns
    -------
    Response
        The response from the Zenodo API after publishing the deposition.
    """
    r = requests.post(
        f"{ZENODO_API_URL}/{deposition_id}/actions/publish",
        params={"access_token": token},
    )
    r.raise_for_status()
    return r


def add_version_row(
    rows, dataset, source, version, tags, url, note, remove_latest=True
):
    """
    Add a new version row to the list of rows.

    Parameters
    ----------
    rows : list of dict
        The existing rows from the versions CSV.
    dataset : str
        The dataset name.
    source : str
        The source of the dataset, likely "archive".
    version : str
        The version name.
    tags : list of str
        The tags for the new version.
    url : str
        The Zenodo URL.
    note : str
        A note for the new version.
    remove_latest : bool, optional
        Whether to remove the 'latest' tag from the existing latest version (default is True).

    Returns
    -------
    list of dict
        The updated rows.
    """

    # Find the index of the row of the latest version for the dataset
    latest_index = next(
        (
            i
            for i, row in enumerate(rows)
            if row["dataset"] == dataset
            and "latest" in row["tags"]
            and row["source"] == source
        ),
        None,
    )

    # If a latest version exists, remove the 'latest' tag from it
    if latest_index is not None and remove_latest:
        rows[latest_index]["tags"] = [
            tag for tag in (rows[latest_index]["tags"]) if tag != "latest"
        ]

    # Add the new version row with 'latest' and 'supported' tags ...
    new_row = {
        "dataset": dataset,
        "source": source,
        "version": version,
        "tags": tags,
        "url": url,
        "note": note,
    }

    # ... either behind the latest version or at the end of the list
    if latest_index is not None:
        # Insert the new row after the latest version
        rows.insert(latest_index + 1, new_row)
    else:
        # Append to the end if no latest version exists
        rows.append(new_row)

    return rows


def get_deposition_by_url(url, token):
    """
    Retrieve a Zenodo deposition by its API URL.

    Parameters
    ----------
    url : str
        The Zenodo deposition API URL.
    token : str
        The Zenodo API access token.

    Returns
    -------
    dict
        The deposition metadata as a dictionary.
    """
    r = requests.get(url, params={"access_token": token})
    r.raise_for_status()
    return r.json()


def zenodo_create_deposition(token, metadata=None):
    """
    Create a new Zenodo deposition.

    Parameters
    ----------
    token : str
        The Zenodo API access token.
    metadata : dict, optional
        The metadata for the deposition.

    Returns
    -------
    dict
        The created deposition metadata.
    """
    headers = {"Content-Type": "application/json"}
    r = requests.post(
        ZENODO_API_URL,
        params={"access_token": token},
        json=metadata or {},
        headers=headers,
    )
    r.raise_for_status()
    return r.json()


def zenodo_update_metadata(deposition_url, token, metadata):
    """
    Update the metadata for a Zenodo deposition.

    Parameters
    ----------
    deposition_url : str
        The Zenodo deposition API URL.
    token : str
        The Zenodo API access token.
    metadata : dict
        The metadata to update.

    Returns
    -------
    dict
        The updated deposition metadata.
    """
    headers = {"Content-Type": "application/json"}
    r = requests.put(
        deposition_url, params={"access_token": token}, json=metadata, headers=headers
    )
    r.raise_for_status()
    return r.json()


def zenodo_publish_deposition(deposition_id, token):
    """
    Publish a Zenodo deposition.

    Parameters
    ----------
    deposition_id : str or int
        The Zenodo deposition ID.
    token : str
        The Zenodo API access token.

    Returns
    -------
    dict
        The published deposition metadata.
    """
    r = requests.post(
        f"{ZENODO_API_URL}/{deposition_id}/actions/publish",
        params={"access_token": token},
    )
    r.raise_for_status()
    return r.json()


def zenodo_new_version(deposition_id, token):
    """
    Create a new version of a Zenodo deposition.

    Parameters
    ----------
    deposition_id : str or int
        The Zenodo deposition ID.
    token : str
        The Zenodo API access token.

    Returns
    -------
    dict
        The new version deposition metadata.
    """
    r = requests.post(
        f"{ZENODO_API_URL}/{deposition_id}/actions/newversion",
        params={"access_token": token},
    )
    r.raise_for_status()
    return r.json()


def extract_zenodo_deposition_url(record_url: str) -> str | None:
    """
    Given a Zenodo record or record file URL, return the corresponding deposition API URL.

    Parameters
    ----------
    record_url : str
        The Zenodo record or record file URL.

    Returns
    -------
    str or None
        The Zenodo deposition API URL, or None if extraction fails.
    """
    # Match both sandbox and production, with or without /files/...
    m = re.match(
        r"https://(sandbox\.)?zenodo\.org/records/(\d+)(/files/.*)?", record_url
    )
    if not m:
        return None
    sandbox, record_id = m.group(1), m.group(2)
    if sandbox:
        return f"https://sandbox.zenodo.org/api/deposit/depositions/{record_id}"
    else:
        return f"https://zenodo.org/api/deposit/depositions/{record_id}"


@app.command()
def main():
    """
    Guide the user through creating a new Zenodo deposition or version.

    This function provides an interactive CLI for creating or releasing new versions
    of datasets on Zenodo, uploading files, and updating the local versions CSV.
    """
    typer.secho("=" * 60, fg=typer.colors.CYAN)
    typer.echo(
        "Welcome! This script helps you create/upload a new dataset version on Zenodo."
    )
    typer.echo(
        "To create/upload a new dataset version, place the dataset in the data/<dataset_name>/archive/<version_name> folder."
    )
    typer.echo(
        "Then run this script; it will guide you through the process and update data/versions.csv accordingly."
    )
    typer.secho("=" * 60, fg=typer.colors.CYAN)

    action = prompt_choice(
        ["release", "create"],
        "Select whether you want to (release) a new version to an existing dataset or (create) a new dataset?",
    )
    token = get_access_token()
    rows = read_versions_csv()

    # All datasets that are recorded in the data/versions.csv and have a 'latest' tag
    latest = get_latest_versions(rows, "archive")

    # Logic is slightly different for 'release' and 'create':
    if action == "release":
        if not latest:
            typer.secho(
                "No datasets with source 'archive' and a Zenodo URL found.",
                fg=typer.colors.RED,
            )
            raise typer.Exit()

        dataset_names = [row["dataset"] for row in latest]
        dataset_name = prompt_choice(
            dataset_names, "Select a dataset to create a new version release for"
        )
    elif action == "create":
        dataset_names = get_potential_datasets()
        dataset_names = [
            ds for ds in dataset_names if ds not in {l["dataset"] for l in latest}
        ]
        if not dataset_names:
            typer.secho(
                "No new datasets found in 'data' that contain an 'archive' subfolder. "
                "Please create a dataset folder with an 'archive' subfolder first in 'data/<dataset_name>/archive/<version_name>'.",
                fg=typer.colors.RED,
            )
            raise typer.Exit()
        dataset_name = prompt_choice(dataset_names, "Select a new dataset to upload")

    # Get potential folders for new version
    typer.echo("Fetching new version folders...")
    potential_folders = get_archive_folders(dataset_name)

    if not potential_folders:
        typer.secho(
            f"No archive folders found in data/{dataset_name}/archive. Please create a folder with the new version files first.",
            fg=typer.colors.RED,
        )
        raise typer.Exit()

    # Exclude all folders for which version entres already exist
    known_versions = get_dataset_versions(rows, dataset_name)
    potential_folders = [
        folder for folder in potential_folders if folder not in known_versions
    ]

    if not potential_folders:
        typer.secho(
            f"All folder in 'data/{dataset_name}/archive already have versions in data/versions.csv. "
            f"To create a new version, please create a new folder with a new version name that includes the new files.",
            fg=typer.colors.RED,
        )
        raise typer.Exit()

    # Ask user which files from which folder to upload and how to call the new version
    folder_to_upload = prompt_choice(
        potential_folders, "Select a folder to upload as new version"
    )
    version_name = typer.prompt(
        "Specify name for new version", default=folder_to_upload
    )

    if action == "release":
        # Retrieve information from previous version for the datset from Zenodo
        row = next(r for r in latest if r["dataset"] == dataset_name)
        deposition_url = extract_zenodo_deposition_url(row["url"])

        if not deposition_url:
            typer.secho(
                f"Could not extract Zenodo deposition API URL from:\n{row['url']}\n"
                "Please check the URL format in your versions.csv.",
                fg=typer.colors.RED,
            )
            raise typer.Exit()

        typer.echo(
            f"Retrieving existing deposition metadata from Zenodo for dataset '{dataset_name}' at {deposition_url}."
        )
        existing_deposition = get_deposition_by_url(deposition_url, token)

        metadata = existing_deposition.get("metadata", {})
        if not metadata:
            typer.secho(
                "No metadata found in the existing deposition. Please check the Zenodo record.",
                fg=typer.colors.RED,
            )
            raise typer.Exit()

        new_metadata = {
            "title": metadata["title"],
            "creators": metadata["creators"],
            "version": version_name,
            "publication_date": datetime.now().strftime("%Y-%m-%d"),
            "description": metadata["description"],
            "license": metadata["license"],
            "upload_type": metadata["upload_type"],
            "access_right": metadata["access_right"],
            "related_identifiers": [
                {
                    "scheme": "doi",
                    "identifier": metadata["doi"],
                    "relation": "isNewVersionOf",
                    "resource_type": metadata["upload_type"],
                }
            ],
        }
    elif action == "create":
        # Create new metadata based on user input
        typer.echo(
            "Creating a new deposition requires metadata input.\n"
            "Please provide the following information for the new dataset."
            "You can also leave fields empty for now and fill them in later on Zenodo before publishing the dataset."
        )
        title = typer.prompt("Title of the new dataset")
        description = typer.prompt("Description of the new dataset")
        license = typer.prompt(
            "License for the new dataset (e.g. CC-BY-4.0)", default="cc-by-4.0"
        )
        upload_type = typer.prompt("Upload type for the new dataset", default="dataset")
        access_right = typer.prompt("Access right for the new dataset", default="open")
        # Prompt for creators and their affiliation (a list of dictionaries with keys "name" and "affiliation") until the user is done
        typer.echo(
            "Please provide the creators of the new dataset.\n"
            "You need to provide at least on creator. You can add multiple creators.\n\n"
            "Leave the name blank when you are done to continue."
        )
        creators = []
        while True:
            name = typer.prompt("Creator name (leave blank to finish)", default="")
            if not name and len(creators) > 0:
                break
            affiliation = typer.prompt("Affiliation (or leave empty)", default="")

            if name:
                creators.append({"name": name, "affiliation": affiliation})

        # Create the new metadata dictionary
        new_metadata = {
            "title": title,
            "description": description,
            "version": version_name,
            "publication_date": datetime.now().strftime("%Y-%m-%d"),
            "license": license,
            "upload_type": upload_type,
            "access_right": access_right,
            "creators": creators,
        }

        # Debug statement:
        typer.echo(f"New metadata created: {new_metadata}")

    # Create deposition and upload files
    typer.echo(
        "Creating new Zenodo deposition with metadata:\n"
        + json.dumps(new_metadata, indent=2)
    )
    deposition = create_zenodo_deposition(
        token,
        ZENODO_API_URL,
        # Provide new metadata baed on the parent deposition
        metadata=new_metadata,
        # all files in the folder
        files=[
            f
            for f in (Path(f"data/{dataset_name}/archive/{folder_to_upload}").glob("*"))
            if f.is_file()
        ],
    )

    # Debug statement:
    typer.echo(
        f"Created new deposition ready for publishing.\n"
        f"You can now review it and make changes at: {deposition.json()['links']['html']}"
    )

    # Confirm publishing
    confirm = typer.confirm("Are you ready to publish the dataset now?", default=True)
    if confirm:
        published_deposition = publish_zenodo_deposition(token, deposition.json()["id"])
        published_deposition.raise_for_status()
        typer.secho(
            f"✅ Dataset '{dataset_name}' with version {version_name} has been successfully published.\n"
            f"➡️ Link to new dataset: {published_deposition.json()['links']['html']}.",
        )
    else:
        published_deposition = None
        typer.secho(
            f"❌ New dataset version not published. "
            f"️️➡️ You can publish it later using the Zenodo web interface at: {deposition.json()['links']['html']}",
        )

    # Update the versions CSV file
    typer.confirm(
        "Proceed with updating the 'data/versions.csv' file with the new version information?",
        default=True,
        abort=True,
    )
    # Debug statement:
    typer.echo(json.dumps(deposition.json(), indent=2))

    previous_version = [row for row in latest if row["dataset"] == dataset_name]
    previous_version = previous_version[0] if previous_version else {}

    if published_deposition:
        tags = ["latest", "supported"]
    elif not published_deposition:
        tags = ["draft", "supported"]

    if published_deposition:
        link = published_deposition.json()["links"]["html"]
    else:
        link = deposition.json()["links"]["html"]

    if previous_version:
        note = previous_version["note"]
    else:
        note = typer.prompt(
            "Please provide a note for the new version in 'data/versions.csv' (or leave empty)",
            default="",
        )

    rows = add_version_row(
        rows,
        dataset_name,
        "archive",
        version_name,
        tags=tags,
        url=link,
        note=note,
        remove_latest=True if published_deposition else False,
    )
    write_versions_csv(rows)
    typer.secho(
        f"✅ 'data/versions.csv' has been updated with: '{dataset_name}', 'archive', '{version_name}', '{tags}, '{link}', '{note}'",
    )

    typer.echo("🎉🎉🎉 All done. Thank you and visit us again!")
    typer.Exit()


if __name__ == "__main__":
    """
    Run the Typer CLI application.
    """
    app()
