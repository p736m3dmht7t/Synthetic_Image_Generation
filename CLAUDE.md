# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview
This is a synthetic image generation project using Python. The project uses a virtual environment located at `.venv/` with Python 3.13.5.  

## Development Environment Setup
- Virtual environment: `.venv/` (Python 3.13.5)
- Activate virtual environment: `source .venv/bin/activate`
- Deactivate virtual environment: `deactivate`

## Common Commands
Since this is a new project, common Python development commands would typically include:
- Install dependencies: `pip install -r requirements.txt` (once requirements.txt exists)
- Run main script: `python main.py` (once main.py exists)
- Run tests: `python -m pytest` (once tests are added)
- Install in development mode: `pip install -e .` (once setup.py/pyproject.toml exists)

## Project Structure
This is a new project directory. The typical structure for a synthetic image generation project would include:
- Source code in a main package directory
- Requirements file for dependencies
- Configuration files for model parameters
- Scripts for training/inference
- Output directories for generated images

## Notes
- Project is currently in initial setup phase
- Virtual environment is configured and ready for development
- Consider adding common ML/AI dependencies like torch, PIL, numpy, etc.

## Important Rules
1. First think through the problem, read the codebase for relevant files, and write a plan to tasks/todo.md.
2. The plan should have a list of todo items that you can check off as you complete them.
3. Before you begin working, check in with me and I will verify the plan.
4. Then, begin working on the todo items, marking them as complete as you go.
5. Please every step of the way just give me a high level explanation of what changes you made
6. Make every task and code change you do as simple as possible. We want to avoid making any massive or complex changes. Every change should impact as little code as possible. Everything is about simplicity.
7. Finally, add a review section to the [todo.md](http://todo.md/) file with a summary of the changes you made and any other relevant information.
8. Follow my preferred commenting style in ./comment_style.md
9. When working on github.com issues, allow the human to test the fix before pushing it.