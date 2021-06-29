
deploy_build <- function(
  pkg = ".",
  branch = "gh-builds",
  remote = "origin"
  ) {
  
  # Create destination directory for built files
  build_dir <- fs::dir_create(fs::file_temp())
  on.exit(fs::dir_delete(build_dir))
  
  # Create desination directory for working tree
  dest_dir <- fs::dir_create(fs::file_temp())
  on.exit(fs::dir_delete(dest_dir))
  
  if (!pkgdown:::git_has_remote_branch(remote, branch)) {
    old_branch <- pkgdown:::git_current_branch()
    
    # If no remote branch, we need to create it
    pkgdown:::git("checkout", "--orphan", branch)
    pkgdown:::git("rm", "-rf", "--quiet", ".")
    pkgdown:::git("commit", "--allow-empty", "-m", sprintf("Initializing %s branch", branch))
    pkgdown:::git("push", remote, paste0("HEAD:", branch))
    
    # checkout the previous branch
    pkgdown:::git("checkout", old_branch)
  }
  
  # Explicitly set the branches tracked by the origin remote.
  # Needed if we are using a shallow clone, such as on travis-CI
  pkgdown:::git("remote", "set-branches", remote, branch)
  
  # Fetch the remote branch
  pkgdown:::git("fetch", remote, branch)
  
  # Create a worktree for the built packages
  pkgdown:::github_worktree_add(dest_dir, remote, branch)
  on.exit(pkgdown:::github_worktree_remove(dest_dir), add = TRUE)
  
  # Get package details
  desc <- pkgdown:::read_desc(pkg)
  package <- desc$get("Package")[[1]]
  version <- desc$get_field("Version")
  
  # Get session details
  R_version <- sessionInfo()$R.version$major
  platform <- sessionInfo()$R.version$platform
  
  # Create a folder for the version number
  R_version_dir <- sprintf("R_v%s", R_version)
  package_filename <- sprintf("%s_%s_%s", package, version, platform)
  package_output_dir <- file.path(dest_dir, R_version_dir)
  dir.create(package_output_dir, showWarnings = FALSE)
  
  # Build the package in the build dir
  pkgbuild::build(
    dest_path = build_dir,
    binary = TRUE
  )
  
  # Move the built package to the output directory
  built_pkg_path <- list.files(build_dir, full.names = T)[1]
  built_pkg_file <- list.files(build_dir)[1]
  built_pkg_name <- gsub(".*\\.", paste0(package_filename, "."), built_pkg_file)
  file.rename(
    built_pkg_path,
    file.path(package_output_dir, built_pkg_name)
  )
  
  # Push the contents of the directory
  commit_message <- paste("build", built_pkg_name)
  pkgdown:::github_push(dest_dir, commit_message, remote, branch)
  invisible()
  
  
}
