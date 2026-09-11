from pathlib import Path

import yaml
from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.statemachine import StringList
from sphinx.util.nodes import nested_parse_with_titles


def _type_name(value):
    return type(value).__name__


class YamlConfigDirective(Directive):
    required_arguments = 1
    has_content = False


    def run(self):
        env = self.state.document.settings.env
        filename = Path(env.srcdir) / self.arguments[0]

        # Tell Sphinx that this document depends on the YAML file. This ensures
        # the page is rebuilt if the YAML changes.
        env.note_dependency(str(filename))

        with filename.open(encoding="utf-8") as f:
            config = yaml.safe_load(f) or {}

        lines = []
        self._render(config, lines)

        container = nodes.section()
        container.document = self.state.document

        content = StringList(
            lines,
            source=str(filename),
        )

        nested_parse_with_titles(
            self.state,
            content,
            container,
        )

        return container.children

    def _render(self, data, lines, prefix=""):
        for key, value in data.items():
            name = f"{prefix}.{key}" if prefix else key

            if isinstance(value, dict):
                self._render(value, lines, name)
                continue

            if isinstance(value, list):
                if not value:
                    lines.extend(
                        [
                            f".. confval:: {name}",
                            "",
                            "   :Type: list",
                            "   :Default: ``[]``",
                            "",
                        ]
                    )
                    continue

                default = value[0]
                alternatives = value[1:]

                lines.extend(
                    [
                        f".. confval:: {name}",
                        "",
                        f"   :Type: {_type_name(default)}",
                        f"   :Default: ``{default!r}``",
                    ]
                )

                if alternatives:
                    lines.extend(
                        [
                            "   :Alternatives:",
                            "",
                        ]
                    )
                    lines.extend(
                        f"      * ``{alternative!r}``" for alternative in alternatives
                    )

                lines.append("")
                continue

            lines.extend(
                [
                    f".. confval:: {name}",
                    "",
                    f"   :Type: {_type_name(value)}",
                    f"   :Default: ``{value!r}``",
                    "",
                ]
            )


def setup(app):
    app.add_directive("yaml-config", YamlConfigDirective)

    return {
        "version": "1.0",
        "parallel_read_safe": True,
    }
