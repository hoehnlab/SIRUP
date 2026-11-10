#' Install the SIRUP skill for an LLM agent
#'
#' `installSIRUPSkill` copies the SIRUP skill into the skills folder of an LLM
#' agent. The skill teaches the agent how to use SIRUP: how to build a SIRUP
#' object, which function answers which question, and which checks to run.
#'
#' @details
#' A skill is a folder with a file named `SKILL.md`. The file starts with a name
#' and a description. At the start of a session, the agent reads the name and the
#' description of each installed skill. It reads the rest of the folder when your
#' request matches a description.
#'
#' SIRUP ships its skill in `inst/skills/sirup`. R copies that folder into the
#' SIRUP folder of your R library when it installs the package. An agent does not
#' read your R library. This function copies the skill to the place where the
#' agent looks.
#'
#' With `install_dir = NULL`, the function uses the skills folder that the agent
#' reads for all your projects:
#'
#' \tabular{ll}{
#'   `agent` \tab Folder for all your projects \cr
#'   `"claude"` \tab `~/.claude/skills` \cr
#'   `"codex"` \tab `~/.agents/skills` \cr
#'   `"antigravity"` \tab `~/.gemini/config/skills` \cr
#' }
#'
#' On Windows, `~` is your user folder (`USERPROFILE`). These folders come from
#' the documentation of each agent. They can change, so check the documentation
#' of your agent if the skill does not appear.
#'
#' To install the skill for one project only, give the skills folder of the
#' project to `install_dir`. For Claude Code this is `.claude/skills` in the
#' project. For Codex and Antigravity it is `.agents/skills`. A relative path
#' starts in the current working directory.
#'
#' The function copies the skill to a folder named `sirup` inside `install_dir`.
#' If that folder exists and holds a `SKILL.md` file, the function removes it
#' first, so that no old files stay. If the folder exists but holds no `SKILL.md`
#' file, the function stops. It does not remove a folder that is not a skill.
#'
#' After the installation, start a new session of the agent and ask which skills
#' it has. Look for `sirup`.
#'
#' @param agent The agent. Use `"claude"` (Claude Code), `"codex"` (Codex), or
#'   `"antigravity"` (Antigravity).
#' @param install_dir The skills folder that receives the skill (default `NULL`).
#'   If `NULL`, the function uses the folder of the agent for all your projects.
#'   Give a path to install the skill for one project, or to use another folder.
#'   The function makes the folder if it does not exist.
#' @param overwrite Replace a skill that is already in `install_dir` (default
#'   `TRUE`)? If `FALSE`, the function stops when the skill exists.
#' @param quiet Hide the message that says where the skill went (default
#'   `FALSE`)?
#'
#' @return The path of the installed skill folder, invisibly.
#'
#' @examples
#' # Install into a temporary folder, so that this example changes nothing on
#' # your computer. To install for real, leave out `install_dir`.
#' skills <- file.path(tempdir(), "agent-skills")
#' installSIRUPSkill("claude", install_dir = skills, quiet = TRUE)
#' list.files(file.path(skills, "sirup"))
#'
#' \dontrun{
#' # Install for all your projects, in the default folder of the agent.
#' installSIRUPSkill("claude")
#' installSIRUPSkill("codex")
#' installSIRUPSkill("antigravity")
#'
#' # Install for one project only.
#' installSIRUPSkill("claude", install_dir = ".claude/skills")
#' installSIRUPSkill("codex", install_dir = ".agents/skills")
#' }
#' @export
installSIRUPSkill <- function(agent = c("claude", "codex", "antigravity"),
                              install_dir = NULL,
                              overwrite = TRUE,
                              quiet = FALSE) {
  agent <- match.arg(agent)
  if (!is.null(install_dir) &&
      (!is.character(install_dir) || length(install_dir) != 1L || is.na(install_dir) ||
       !nzchar(install_dir))) {
    rlang::abort("`install_dir` must be NULL or a single non-empty folder path.")
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    rlang::abort("`overwrite` must be TRUE or FALSE.")
  }

  skills_dir <- if (is.null(install_dir)) {
    .sirup_skill_default_dir(agent)
  } else {
    path.expand(install_dir)
  }

  source <- system.file("skills", "sirup", package = "SIRUP")
  if (!nzchar(source) || !file.exists(file.path(source, "SKILL.md"))) {
    rlang::abort(c(
      "This installation of SIRUP has no skill folder.",
      "i" = "Install a SIRUP version that ships `inst/skills/sirup`, then run this function again."
    ))
  }

  target <- file.path(skills_dir, "sirup")
  if (file.exists(target)) {
    if (!isTRUE(overwrite)) {
      rlang::abort(c(
        paste0("The skill folder already exists: ", target),
        "i" = "Use `overwrite = TRUE` to replace it."
      ))
    }
    if (!file.exists(file.path(target, "SKILL.md"))) {
      rlang::abort(c(
        paste0("A folder named `sirup` exists in `install_dir` and it is not a skill: ", target),
        "i" = "The function does not remove it. Move it, or choose another `install_dir`."
      ))
    }
    unlink(target, recursive = TRUE)
  }

  dir.create(skills_dir, recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(source, skills_dir, recursive = TRUE)
  if (!isTRUE(copied) || !file.exists(file.path(target, "SKILL.md"))) {
    rlang::abort(paste0("The function could not copy the skill to: ", skills_dir))
  }

  if (!isTRUE(quiet)) {
    cli::cli_inform(c(
      "v" = "Installed the SIRUP skill for {.val {agent}} in {.path {target}}.",
      "i" = .sirup_skill_next_step(agent)
    ))
  }
  invisible(target)
}

# The user's home folder. On Windows, R's `~` is the Documents folder, so use USERPROFILE.
.sirup_home_dir <- function() {
  home <- Sys.getenv(if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME")
  if (!nzchar(home)) home <- path.expand("~")
  home
}

# The skills folder that an agent reads for all projects.
.sirup_skill_default_dir <- function(agent, home = .sirup_home_dir()) {
  switch(agent,
    claude = file.path(home, ".claude", "skills"),
    codex = file.path(home, ".agents", "skills"),
    antigravity = file.path(home, ".gemini", "config", "skills"),
    rlang::abort(paste0("Unknown agent: ", agent))
  )
}

.sirup_skill_next_step <- function(agent) {
  switch(agent,
    claude = "Start a new Claude Code session and ask which skills it has. Look for `sirup`.",
    codex = "Restart Codex if the skill does not appear. Run `/skills` to list the skills.",
    antigravity = paste0(
      "Start a new session and ask which skills the agent has. ",
      "In the Antigravity command-line tool, run `/skills reload`."
    )
  )
}
