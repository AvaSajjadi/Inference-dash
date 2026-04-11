#!/bin/bash
set -e

echo "Installing Heroku CLI..."
curl https://cli-assets.heroku.com/install.sh | sh

echo "Verifying installation..."
heroku --version

echo "Logging in to Heroku..."
heroku login

echo "Creating Heroku app..."
heroku create

echo "Done! Share your app name with Claude."
