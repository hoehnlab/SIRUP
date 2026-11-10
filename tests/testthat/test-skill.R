skip_if_not(nzchar(system.file("skills", "sirup", package = "SIRUP")),
            "This installation of SIRUP has no skill folder.")

test_that("each agent has its own default skills folder", {
  expect_equal(.sirup_skill_default_dir("claude", home = "/h"), file.path("/h", ".claude", "skills"))
  expect_equal(.sirup_skill_default_dir("codex", home = "/h"), file.path("/h", ".agents", "skills"))
  expect_equal(
    .sirup_skill_default_dir("antigravity", home = "/h"),
    file.path("/h", ".gemini", "config", "skills")
  )
  expect_error(.sirup_skill_default_dir("unknown", home = "/h"), "Unknown agent")
})

test_that("installSIRUPSkill copies the skill into install_dir for every agent", {
  for (agent in c("claude", "codex", "antigravity")) {
    skills <- withr::local_tempdir()
    path <- installSIRUPSkill(agent, install_dir = skills, quiet = TRUE)
    expect_equal(normalizePath(path), normalizePath(file.path(skills, "sirup")))
    expect_true(file.exists(file.path(path, "SKILL.md")))
    expect_true(dir.exists(file.path(path, "references")))
  }
})

test_that("installSIRUPSkill makes an install_dir that does not exist", {
  skills <- file.path(withr::local_tempdir(), "a", "b", "skills")
  installSIRUPSkill("claude", install_dir = skills, quiet = TRUE)
  expect_true(file.exists(file.path(skills, "sirup", "SKILL.md")))
})

test_that("the default folder follows the home folder of the user", {
  home <- withr::local_tempdir()
  withr::local_envvar(HOME = home, USERPROFILE = home)
  installSIRUPSkill("codex", quiet = TRUE)
  expect_true(file.exists(file.path(home, ".agents", "skills", "sirup", "SKILL.md")))
  installSIRUPSkill("antigravity", quiet = TRUE)
  expect_true(file.exists(file.path(home, ".gemini", "config", "skills", "sirup", "SKILL.md")))
  installSIRUPSkill("claude", quiet = TRUE)
  expect_true(file.exists(file.path(home, ".claude", "skills", "sirup", "SKILL.md")))
})

test_that("an update removes the files that the old copy had", {
  skills <- withr::local_tempdir()
  installSIRUPSkill("claude", install_dir = skills, quiet = TRUE)
  stale <- file.path(skills, "sirup", "references", "stale.md")
  writeLines("old", stale)
  installSIRUPSkill("claude", install_dir = skills, quiet = TRUE)
  expect_false(file.exists(stale))
  expect_true(file.exists(file.path(skills, "sirup", "SKILL.md")))
})

test_that("overwrite = FALSE stops when the skill exists", {
  skills <- withr::local_tempdir()
  installSIRUPSkill("claude", install_dir = skills, quiet = TRUE)
  expect_error(
    installSIRUPSkill("claude", install_dir = skills, overwrite = FALSE, quiet = TRUE),
    "already exists"
  )
  expect_true(file.exists(file.path(skills, "sirup", "SKILL.md")))
})

test_that("the function does not remove a folder that is not a skill", {
  skills <- withr::local_tempdir()
  other <- file.path(skills, "sirup")
  dir.create(other)
  writeLines("mine", file.path(other, "notes.txt"))
  expect_error(installSIRUPSkill("claude", install_dir = skills, quiet = TRUE), "not a skill")
  expect_true(file.exists(file.path(other, "notes.txt")))
})

test_that("installSIRUPSkill checks its arguments", {
  expect_error(installSIRUPSkill("gemini"))
  expect_error(installSIRUPSkill("claude", install_dir = c("a", "b")), "install_dir")
  expect_error(installSIRUPSkill("claude", install_dir = ""), "install_dir")
  expect_error(installSIRUPSkill("claude", install_dir = 1), "install_dir")
  expect_error(installSIRUPSkill("claude", install_dir = withr::local_tempdir(), overwrite = NA), "overwrite")
})

test_that("installSIRUPSkill says where the skill went unless quiet is TRUE", {
  skills <- withr::local_tempdir()
  expect_message(installSIRUPSkill("codex", install_dir = skills), "Installed the SIRUP skill")
  expect_no_message(installSIRUPSkill("codex", install_dir = skills, quiet = TRUE))
})
