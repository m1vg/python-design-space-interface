define([
    "nbextensions/widgets/widgets/js/widget",
    "jqueryui",
    "underscore",
    "base/js/utils",
    "bootstrap",
], function(widget, $, _, utils) {
    "use strict";

    var PopupView = widget.DOMWidgetView.extend({
        initialize: function() {
            /**
             * Public constructor
             */
            PopupView.__super__.initialize.apply(this, arguments);
            this.children_views = new widget.ViewList(this.add_child_model, null, this);
            this.listenTo(this.model, 'change:children', function(model, value) {
                this.children_views.update(value);
            }, this);
            this.listenTo(this.model, 'change:overflow_x', this.update_overflow_x, this);
            this.listenTo(this.model, 'change:overflow_y', this.update_overflow_y, this);
            this.listenTo(this.model, 'change:box_style', this.update_box_style, this);
            this.$el.resizable();
            this.$el.css('background', 'white');
            this.$el.css('z-index', '4');
        },

        update_attr: function(name, value) {
            /**
             * Set a css attr of the widget view.
             */
            this.$box.css(name, value);
        },

        render: function() {
            /**
             * Called when view is rendered.
             */
            this.$box = this.$el;
            this.$el.append($('<div class="draggable" style="height:15px; border-style:groove;"></div>'))
            this.$box.draggable({ handle: ".draggable" });
            this.$box.addClass('widget-box');
            this.children_views.update(this.model.get('children'));
            this.update_overflow_x();
            this.update_overflow_y();
            this.update_box_style('');
        },

        update_overflow_x: function() {
            /**
             * Called when the x-axis overflow setting is changed.
             */
            this.$box.css('overflow-x', this.model.get('overflow_x'));
        },

        update_overflow_y: function() {
            /**
             * Called when the y-axis overflow setting is changed.
             */
            this.$box.css('overflow-y', this.model.get('overflow_y'));
        },

        update_box_style: function(previous_trait_value) {
            var class_map = {
                success: ['alert', 'alert-success'],
                info: ['alert', 'alert-info'],
                warning: ['alert', 'alert-warning'],
                danger: ['alert', 'alert-danger']
            };
            this.update_mapped_classes(class_map, 'box_style', previous_trait_value, this.$box[0]);
        },

        add_child_model: function(model) {
            /**
             * Called when a model is added to the children list.
             */
            var that = this;
            var dummy = $('<div/>');
            that.$box.append(dummy);
            return this.create_child_view(model).then(function(view) {
                dummy.replaceWith(view.el);
                // Trigger the displayed event of the child view.
                that.displayed.then(function() {
                    view.trigger('displayed');
                });
                return view;
            }).catch(utils.reject("Couldn't add child view to box", true));
        },

        remove: function() {
            /**
             * We remove this widget before removing the children as an optimization
             * we want to remove the entire container from the DOM first before
             * removing each individual child separately.
             */
            PopupView.__super__.remove.apply(this, arguments);
            this.children_views.remove();
        },
    });



    return {
        PopupView: PopupView
    }
});