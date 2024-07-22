#!/bin/sh

cp --recursive ../pypsa-eur/* ./ 
echo 'solved tutorial files are copied to workspace'

cp .devcontainer/welcome-message.txt /usr/local/etc/vscode-dev-containers/first-run-notice.txt