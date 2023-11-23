# Welcome to lavaan contributing guide <!-- omit in toc -->

Thank you for investing your time in contributing to the lavaan project!

Read our [Code of Conduct](./code_of_conduct.md) to keep our community approachable and respectable.

In this small guide you will get an overview of the contribution workflow from
opening an issue, creating a pull request (=PR), reviewing, and merging the PR.

## New contributor guide

Here are some resources to help you get started with open source contributions:

- [Finding ways to contribute to open source on GitHub](https://docs.github.com/en/get-started/exploring-projects-on-github/finding-ways-to-contribute-to-open-source-on-github)
- [Set up Git](https://docs.github.com/en/get-started/quickstart/set-up-git)
- [GitHub flow](https://docs.github.com/en/get-started/quickstart/github-flow)
- [Collaborating with pull requests](https://docs.github.com/en/github/collaborating-with-pull-requests)

## Getting started

To understand how lavaan works and how the code base is organized, see this
[understanding_lavaan_internals.R](https://github.com/yrosseel/lavaan/blob/master/inst/understanding_lavaan_internals.R) document. When in doubt, feel
free to contact the main developer.

### Issues

#### Create a new issue

If you have discovered a potential issue, [check if the issue already exists](https://github.com/yrosseel/lavaan/issues). If not, you can open a new issue.

If you're unable to find an open issue addressing the problem, [open an new
issue](https://github.com/yrosseel/lavaan/issues/new). Be sure to include a
title and clear description, as much relevant information as possible, and a
code sample or an executable test case demonstrating the expected behavior that
is not occurring, or the feature that you are requesting.

#### Solve an issue

Scan through our [existing issues](https://github.com/yrosseel/lavaan/issues)
to find one that interests you. As a general rule, we don’t assign issues to
anyone. If you find an issue to work on, you are welcome to open a PR with a
fix.

### Make Changes

1. Fork the repository.
- Using GitHub Desktop:
  - [Getting started with GitHub Desktop](https://docs.github.com/en/desktop/installing-and-configuring-github-desktop/getting-started-with-github-desktop) will guide you through setting up Desktop.
  - Once Desktop is set up, you can use it to [fork the repo](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/cloning-and-forking-repositories-from-github-desktop)!

- Using the command line:
  - [Fork the repo](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo#fork-an-example-repository) so that you can make your changes without affecting the original project until you're ready to merge them.

2. Create a working branch and start with your changes!

### Commit your update

Commit the changes once you are happy with them.

### Pull Request

When you're finished with the changes, create a pull request, also known as a PR.
- Fill the "Ready for review" template so that we can review your PR. This template helps reviewers understand your changes as well as the purpose of your pull request.
- Don't forget to [link PR to issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue) if you are solving one.
- Enable the checkbox to [allow maintainer edits](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/allowing-changes-to-a-pull-request-branch-created-from-a-fork) so the branch can be updated for a merge.
Once you submit your PR, the lavaan maintainer will review your proposal. We may ask questions or request additional information.
- We may ask for changes to be made before a PR can be merged.
- As you update your PR and apply changes, mark each conversation as [resolved](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/commenting-on-a-pull-request#resolving-conversations).
- If you run into any merge issues, checkout this [git tutorial](https://github.com/skills/resolve-merge-conflicts) to help you resolve merge conflicts and other issues.

