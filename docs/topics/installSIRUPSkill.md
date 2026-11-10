# installSIRUPSkill

## Install the SIRUP skill for an LLM agent

## Description

`installSIRUPSkill` copies the SIRUP skill into the skills folder of an LLM agent. The skill teaches the agent how to use SIRUP: how to build a SIRUP object, which function answers which question, and which checks to run.

## Usage

``` R
installSIRUPSkill(
  agent = c("claude", "codex", "antigravity"),
  install_dir = NULL,
  overwrite = TRUE,
  quiet = FALSE
)
```

## Arguments

`agent`

:   The agent. Use `"claude"` (Claude Code), `"codex"` (Codex), or `"antigravity"` (Antigravity).

`install_dir`

:   The skills folder that receives the skill (default `NULL`). If `NULL`, the function uses the folder of the agent for all your projects. Give a path to install the skill for one project, or to use another folder. The function makes the folder if it does not exist.

`overwrite`

:   Replace a skill that is already in `install_dir` (default `TRUE`)? If `FALSE`, the function stops when the skill exists.

`quiet`

:   Hide the message that says where the skill went (default `FALSE`)?

## Details

A skill is a folder with a file named `SKILL.md`. The file starts with a name and a description. At the start of a session, the agent reads the name and the description of each installed skill. It reads the rest of the folder when your request matches a description.

SIRUP ships its skill in `inst/skills/sirup`. R copies that folder into the SIRUP folder of your R library when it installs the package. An agent does not read your R library. This function copies the skill to the place where the agent looks.

With `install_dir = NULL`, the function uses the skills folder that the agent reads for all your projects:

`agent`
:   Folder for all your projects

`"claude"`
:   `⁠~/.claude/skills⁠`

`"codex"`
:   `⁠~/.agents/skills⁠`

`"antigravity"`
:   `⁠~/.gemini/config/skills⁠`

On Windows, `~` is your user folder (`USERPROFILE`). These folders come from the documentation of each agent. They can change, so check the documentation of your agent if the skill does not appear.

To install the skill for one project only, give the skills folder of the project to `install_dir`. For Claude Code this is `.claude/skills` in the project. For Codex and Antigravity it is `.agents/skills`. A relative path starts in the current working directory.

The function copies the skill to a folder named `sirup` inside `install_dir`. If that folder exists and holds a `SKILL.md` file, the function removes it first, so that no old files stay. If the folder exists but holds no `SKILL.md` file, the function stops. It does not remove a folder that is not a skill.

After the installation, start a new session of the agent and ask which skills it has. Look for `sirup`.

## Value

The path of the installed skill folder, invisibly.

## Examples

``` R
# Install into a temporary folder, so that this example changes nothing on
# your computer. To install for real, leave out `install_dir`.
skills <- file.path(tempdir(), "agent-skills")
installSIRUPSkill("claude", install_dir = skills, quiet = TRUE)
list.files(file.path(skills, "sirup"))

## Not run: 
# Install for all your projects, in the default folder of the agent.
installSIRUPSkill("claude")
installSIRUPSkill("codex")
installSIRUPSkill("antigravity")

# Install for one project only.
installSIRUPSkill("claude", install_dir = ".claude/skills")
installSIRUPSkill("codex", install_dir = ".agents/skills")

## End(Not run)
```
