# Publishing the prepared code with GitHub Desktop

These steps apply when the authors decide to publish the scripts and choose a
license. They do not upload anything automatically. A public source archive can
be shared before a full clean-environment rerun, provided its limitations are
described accurately.

1. Install GitHub Desktop and sign in to the account that will own the code.
2. In GitHub Desktop, use **File > Add local repository** and select the
   `11_Code_Release_Staging` directory. If prompted, choose **Create a
   repository here**. Use a short repository name such as
   `gastric-plat-fibroblast-macrophage`.
3. Review the file list. The repository should contain `code/`, `gene_sets/`,
   `resources/`, `docs/`, `README.md`, and `.gitignore`. Read
   `resources/README.md` and keep its MSigDB and NicheNet-v2 attribution with
   the copied resources. Confirm that `_stage_code.py` and
   `manuscript_revisions_draft/` are ignored. Do not add the full research
   project, raw data, local environments, submission documents or signed forms.
4. Commit the files locally, then select **Publish repository**. Choose
   **Public** only when the authors have approved release.
5. Open the GitHub repository URL and verify that README, script manifest,
   resource manifest, primary signature TSV and seed registry are visible.
   Record the commit ID or create a version tag/release. Add that exact URL
   and version to the manuscript Data sharing or Code availability statement.

For a command-line workflow, the same directory can be initialized with Git,
committed, linked to a newly created empty GitHub repository, and pushed.
Never upload the whole `ECM_CFD_CXCL12_BCELL_GC_2026` project directory.
