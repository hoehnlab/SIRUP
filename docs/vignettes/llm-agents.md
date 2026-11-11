Use SIRUP with your LLM agent
================

An LLM agent can write and run SIRUP code for you, and it can write
better code when it knows how SIRUP works. SIRUP provides agents with:
\* a packaged skill \* the llms files, which are text and JSON copies of
the documentation.

## How to teach your agent about SIRUP

| Method | What is it? | When to choose this? |
|----|----|----|
| The SIRUP skill | An agent skill | Your agent supports skills. Claude Code, Codex, and Antigravity (and likely others) do. |
| The llms files | Text and JSON copies of the documentation | Your agent does not support skills or you do not wish to use one. |
| The readthedocs pages | The pages of this site | You don’t want to upload files to your agent, and it can browse the web. |

## Install the skill

Agent Skills is an open format created by Anthropic, and it is descrbied
at [agentskills.io](https://agentskills.io). A skill is a folder with a
file named `SKILL.md` that explains how to accomplish a task in a
certain way to an agent. The file starts with a name and a description.
At the start of a session, the agent reads the name and the description
of each installed skill. It reads the rest of the folder only when you
request that it uses the skill.

SIRUP ships its skill inside the R package, in `inst/skills/sirup`,
which is automatically copied into the package folder when you install a
package in R. This puts the skill in your R library, so we need to put
it in a place that your agent will look. For convenience, the
`installSIRUPSkill()` function copies the skill to the appropriate
location for Claude, Codex, and Antigravity agents.

``` r
library(SIRUP)

installSIRUPSkill(agent = "claude")        # Claude Code
installSIRUPSkill(agent = "codex")         # Codex
installSIRUPSkill(agent = "antigravity")   # Antigravity
```

Run the line for your agent. By default, the function installs the skill
in the folder that the agent reads for all your projects.

| `agent`                       | Skill folder for all projects |     |
|-------------------------------|-------------------------------|-----|
| `"claude"` (Claude Code)      | `~/.claude/skills/`           |     |
| `"codex"` (Codex)             | `~/.agents/skills/`           |     |
| `"antigravity"` (Antigravity) | `~/.gemini/config/skills/`    |     |

If you only wish to install a skill for a particular project, use the
`install_dir` argument. Claude agents read `.claude/skills/` in the
project, and Codex and Antigravity agents read `.agents/skills/` in the
project.

``` r
installSIRUPSkill(agent = "claude", install_dir = ".claude/skills")
installSIRUPSkill(agent = "codex", install_dir = ".agents/skills")
installSIRUPSkill(agent = "antigravity", install_dir = ".agents/skills")
```

<div markdown="1"
style="background: #9e584a; padding: 12px; margin-bottom: 24px; line-height: 24px;">

<p style="margin: -12px -12px 12px; padding: 6px 12px; background: #8b3620; color: #ffffff; font-weight: 700;">

<b>WARNING</b>
</p>

These skill locations are correct as of late 2026. Updated architectures
or other agent platforms may change them.

</div>

You can also copy the skill by hand.
`system.file("skills", "sirup", package = "SIRUP")` returns the path of
the skill folder. Copy that folder into the skills folder of your agent.

### Check that your agent sees the skill

Start a new session. Ask the agent which skills it has, and look for
`sirup`. Alternatively use `/skills` in a command line interface to list
all available skills. In the Antigravity editor, open the Customizations
menu in the agent panel. If you add a skill while an agent runs, restart
the agent or reload skills (as with `/reload-skills` for Claude) first
to see the new skill.

### Update the skill

New SIRUP versions automatically keep the skill updated. After you
update SIRUP, run `installSIRUPSkill()` once again. It replaces the old
copy by default.

## Use the llms files instead

The SIRUP repository holds four files specfically for LLMs.

| File             | Content                                                 |
|------------------|---------------------------------------------------------|
| `llms.txt`       | A short Markdown index of the documentation, with links |
| `llms.json`      | The same short index in JSON                            |
| `llms-full.txt`  | The whole set of documentation in one text file         |
| `llms-full.json` | The whole set of documentation in JSON                  |

Give an agent that has no skill support the address of `llms.json` or
`llms.txt`. Tell it to read the index and then to fetch the page that it
needs. The full files are large, about 400 KB each, so chat windows may
not be able to hold them. We recommend that you start with the short
index.

## Example requests

Use these requests with the skill installed. Replace the file names with
your own.

```
I have a Seurat object in `spatial.rds` and a corresponding spatial AIRR rearrangement in `bcr.tsv.gz`.
Please build a SIRUP object and then build clones. Tell me how many sequences you removed and why.

Which spatial statistic should I use to test whether one clone stays in one place?

Load the region polygons in polygons.csv. Which clones dominate each region?

Build lineage trees for the three largest clones in region 5. Tell me what the run needs
before you start.

Write a Slurm script that runs `processSpatialSample()` on my demultiplexed reads.

This error came from `createSIRUP()`: <paste the error>. Can you please explain it and fix the call?
```

## Work with an agent with care

- **Read the code before you run it.** Ask the agent to explain each
  step.
- **Protect private data.** A hosted agent sends your requests, and the
  files that it reads, to its provider. Check your data agreement before
  you share patient data. The example data that comes with SIRUP is for
  tutorials.
- **Ask for checks.** Ask the agent to print the receptor count before
  and after each change, and to run `validateSIRUP()`.
- **Ask what the agent did not run.** An agent can describe a result
  that it did not produce.
- **Agents cannot draw ROIs.** Unless you give your agent a skill that
  allows it to click, you will have to draw regions with the polygon
  selection RShiny app. These can be subjective at times, so we
  recommend not passing this task to an agent.

**An agent does not replace an expert review of your results.**
